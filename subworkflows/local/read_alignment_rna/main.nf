//
// Align RNA reads and lift alignments back to genome coordinates
//

include { BWAMEM2_ALIGN_RNA } from '../../../modules/local/bwa-mem2/mem/rna/main'
include { FASTP_SPLIT       } from '../../../modules/local/fastp/split/main'
include { TARS              } from '../../../modules/local/tars/main'

workflow READ_ALIGNMENT_RNA {
    take:
    // Sample data
    ch_inputs            // channel: [mandatory] [ meta ]
    ch_fastq             // channel: [mandatory] [ meta, fastq_info, fastq_fwd, fastq_rev ]

    // Reference data
    genome_fasta         // channel: [mandatory] /path/to/genome_fasta
    genome_version       // channel: [mandatory] genome version
    genome_fai           // channel: [mandatory] /path/to/genome_fai
    genome_dict          // channel: [mandatory] /path/to/genome_dict
    genome_bwamem2_index // channel: [mandatory] /path/to/genome_bwa-mem2_index_dir/
    contigs_mapping_rna  // channel: [mandatory] /path/to/contigs_mapping_rna
    unmap_regions_rna    // channel: [mandatory] /path/to/unmap_regions_rna

    // Params
    max_fastq_records    // numeric: [optional]  max number of FASTQ records per split

    main:
    //
    // STEP: Handle inputs
    //
    // Sort inputs
    // runnable: channel: [ meta, fastq_info, fastq_fwd, fastq_rev ]
    // skip: channel: [ meta ]

    // NOTE(LN): RNA is only ever sequenced for the tumor sample
    ch_inputs_rna_sorted = ch_fastq
        .branch { meta, fastq_info, fastq_fwd, fastq_rev ->
            def has_inputs = fastq_fwd && fastq_rev
            runnable: fastq_info.sample_type == 'tumor' && has_inputs
            skip: fastq_info.sample_type == 'tumor' && ! has_inputs
              return meta
        }

    // Create FASTQ input channel
    // channel: [ meta_fastq, fastq_fwd, fastq_rev ]
    ch_fastq_inputs = channel.empty()
        .mix(ch_inputs_rna_sorted.runnable)
        .map { meta, fastq_info, fastq_fwd, fastq_rev ->

            // NOTE(SW): initial map sets defaults and conventional ordering of selected fields, merging then overwrites / adds while preserving order
            def rg_id = [fastq_info.sample_id, fastq_info.library_id, fastq_info.lane, fastq_info.flowcell].findAll().join('.')
            def rg_entries = [ID: rg_id, SM: fastq_info.sample_id, LB: fastq_info.library_id] + fastq_info.rg_fields
            def rg_line = '@RG\\t' + rg_entries.collect { k, v -> "${k}:${v}" }.join('\\t')

            def meta_fastq = [
                key: meta.group_id,
                id: "${meta.group_id}_${fastq_info.sample_id}_${fastq_info.library_id}_${fastq_info.lane}",
                rg_line: rg_line,
                sample_id: fastq_info.sample_id,
                library_id: fastq_info.library_id,
                lane: fastq_info.lane,
                output_file_id: rg_id,
                sample_type: fastq_info.sample_type,
            ]

            if (fastq_info.flowcell) {
                meta_fastq.id = "${meta_fastq.id}_${fastq_info.flowcell}"
            }

            return [meta_fastq, fastq_fwd, fastq_rev]

        }

    //
    // MODULE: fastp
    //
    // Split FASTQ into chunks if requested for distributed processing
    // channel: [ meta_fastq_ready, fastq_fwd, fastq_fwd ]
    ch_fastqs_ready = channel.empty()
    // NOTE(SW): required for strict syntax without params block declaration
    if (max_fastq_records.toInteger() > 0) {

        // Run process
        FASTP_SPLIT(
            ch_fastq_inputs,
            // NOTE(SW): required for strict syntax without params block declaration
            max_fastq_records.toInteger(),
        )

        // NOTE(LN): the transpose operator pairs the R1 and R2 chunks by index, and also covers the single chunk case
        // where fastp emits one file per read rather than a list
        ch_fastqs_ready = FASTP_SPLIT.out[0]
            .transpose()
            .map { meta_fastq, fwd, rev ->

                def split_fwd = fwd.name.replaceAll('\\..+$', '')
                def split_rev = rev.name.replaceAll('\\..+$', '')

                assert split_fwd == split_rev

                // NOTE(SW): split allows meta_fastq_ready to be unique, which is required during reunite below
                def meta_fastq_ready = meta_fastq + [id: "${meta_fastq.id}_${split_fwd}", split: split_fwd]

                return [meta_fastq_ready, fwd, rev]
            }

    } else {

        ch_fastqs_ready = ch_fastq_inputs
            .map { meta_fastq, fastq_fwd, fastq_rev ->

                def meta_fastq_ready = meta_fastq + [split: null]

                return [meta_fastq_ready, fastq_fwd, fastq_rev]
            }

    }

    //
    // MODULE: BWA-MEM2
    //
    // Create process input channel
    // channel: [ meta_bwamem2, fastq_fwd, fastq_rev ]
    ch_bwamem2_inputs = ch_fastqs_ready
        .map { meta_fastq_ready, fastq_fwd, fastq_rev ->
            def meta_bwamem2 = meta_fastq_ready.clone()
            return [meta_bwamem2, fastq_fwd, fastq_rev]
        }

    // Run process
    BWAMEM2_ALIGN_RNA(
        ch_bwamem2_inputs,
        genome_fasta,
        genome_bwamem2_index,
    )

    // Reunite BAMs
    // Count expected BAMs per sample for non-blocking groupTuple op
    // channel: [ meta_group, group_size ]
    ch_sample_fastq_counts = ch_bwamem2_inputs
        .map { meta_bwamem2, _reads_fwd, _reads_rev ->

            def meta_group = [
                key: meta_bwamem2.key,
                sample_type: meta_bwamem2.sample_type,
            ]

            return [meta_group, meta_bwamem2]
        }
        .groupTuple()
        .map { meta_group, metas_bwamem2 -> return [meta_group, metas_bwamem2.size()] }

    // Now, group with expected size
    // NOTE(LN): RNA alignments are name-grouped and therefore unindexed, so no index is carried here
    // channel: [ meta_group, [aln, ...] ]
    ch_alns_united = ch_sample_fastq_counts
        // channel: [ [ meta_group, count ], [ meta_group, aln ] ]
        .cross(
            // First element to match meta_group above for `cross`
            channel.topic('bwamem2_align_rna_bam').map { meta_bwamem2, aln -> [[key: meta_bwamem2.key, sample_type: meta_bwamem2.sample_type], aln] }
        )
        .map { count_tuple, inputs_tuple ->
            def group_size = count_tuple[1]
            def (meta_group, aln) = inputs_tuple

            return tuple(groupKey(meta_group, group_size), aln)
        }
        .groupTuple()

    //
    // MODULE: TARS
    //
    // Create process input channel, restoring original meta to source the RNA sample name
    // channel: [ meta_tars, [aln, ...] ]
    ch_tars_inputs = WorkflowOncoanalyser.restoreMeta(ch_alns_united, ch_inputs)
        .map { meta, alns ->

            def sample_id = Utils.getTumorRnaSampleName(meta)

            def meta_tars = [
                key: meta.group_id,
                id: "${meta.group_id}_${sample_id}",
                sample_id: sample_id,
            ]

            return [meta_tars, alns]
        }

    // Run process
    TARS(
        ch_tars_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        contigs_mapping_rna,
        unmap_regions_rna,
    )

    //
    // STEP: Handle outputs
    //
    // Set outputs, restoring original meta
    // NOTE(LN): Tars emits one BAM per sample, but it is carried as a list so that the output shape matches the DNA
    // read alignment subworkflow
    // channel: [ meta, [aln, ...], [idx, ...] ]
    ch_outputs_rna = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('tars_bam'), ch_inputs)
                .map { meta, aln, idx -> [meta, [aln], [idx]] },
            ch_inputs_rna_sorted.skip.unique().map { meta -> [meta, [], []] },
        )

    emit:
    rna = ch_outputs_rna // channel: [ meta, [aln, ...], [idx, ...] ]
}
