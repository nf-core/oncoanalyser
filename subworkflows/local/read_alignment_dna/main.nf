//
// Align DNA reads
//

include { BWAMEM2_ALIGN } from '../../../modules/local/bwa-mem2/mem/main'
include { FASTP_SPLIT   } from '../../../modules/local/fastp/split/main'

workflow READ_ALIGNMENT_DNA {
    take:
    // Sample data
    ch_inputs            // channel: [mandatory] [ meta ]
    ch_fastq             // channel: [mandatory] [ meta, fastq_info, fastq_fwd, fastq_rev ]

    // Reference data
    genome_fasta         // channel: [mandatory] /path/to/genome_fasta
    genome_bwamem2_index // channel: [mandatory] /path/to/genome_bwa-mem2_index_dir/

    // Params
    max_fastq_records    // numeric: [optional]  max number of FASTQ records per split

    main:
    //
    // STEP: Handle inputs
    //
    // Sort inputs
    // runnable: channel: [ meta, fastq_info, fastq_fwd, fastq_rev ]
    // skip: channel: [ meta ]

    // NOTE(SW): the logic block diverges from the standard to accommodate the nature of inputs
    //   - existing input inferred from upstream check
    //   - channel elements are each a single FASTQ pair
    //   - hence case meta maps to pairs as one-to-many
    //   - skip conditional is adjusted for this albiet emits duplicate case meta
    //   - case meta in skip channel deduplicated at join below

    ch_inputs_tumor_sorted = ch_fastq
        .branch { meta, fastq_info, fastq_fwd, fastq_rev ->
            def has_inputs = fastq_fwd && fastq_rev
            runnable: fastq_info.sample_type == 'tumor' && has_inputs
            skip: ! Utils.hasTumorDnaFastq(meta)
              return meta
        }

    ch_inputs_normal_sorted = ch_fastq
        .branch { meta, fastq_info, fastq_fwd, fastq_rev ->
            def has_inputs = fastq_fwd && fastq_rev
            runnable: fastq_info.sample_type == 'normal' && has_inputs
            skip: ! Utils.hasNormalDnaFastq(meta)
              return meta
        }

    ch_inputs_donor_sorted = ch_fastq
        .branch { meta, fastq_info, fastq_fwd, fastq_rev ->
            def has_inputs = fastq_fwd && fastq_rev
            runnable: fastq_info.sample_type == 'donor' && has_inputs
            skip: ! Utils.hasDonorDnaFastq(meta)
              return meta
        }

    // Create FASTQ input channel
    // channel: [ meta_fastq, fastq_fwd, fastq_rev ]
    ch_fastq_inputs = channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable,
            ch_inputs_normal_sorted.runnable,
            ch_inputs_donor_sorted.runnable,
        )
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

        // Now prepare according to FASTQs splitting
        ch_fastqs_ready = channel.topic('fastp_split_fastq')
            .flatMap { meta_fastq, reads_fwd_input, reads_rev_input ->

                def reads_fwd = reads_fwd_input instanceof List ? reads_fwd_input : [reads_fwd_input]
                def reads_rev = reads_rev_input instanceof List ? reads_rev_input : [reads_rev_input]

                def data = [reads_fwd, reads_rev]
                    .transpose()
                    .collect { fwd, rev ->

                        def split_fwd = fwd.name.replaceAll('\\..+$', '')
                        def split_rev = rev.name.replaceAll('\\..+$', '')

                        assert split_fwd == split_rev

                        // NOTE(SW): split allows meta_fastq_ready to be unique, which is required during reunite below
                        def meta_fastq_ready = meta_fastq + [id: "${meta_fastq.id}_${split_fwd}", split: split_fwd]

                        return [meta_fastq_ready, fwd, rev]
                    }

                return data
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
    BWAMEM2_ALIGN(
        ch_bwamem2_inputs,
        genome_fasta,
        genome_bwamem2_index,
    )

    // Reunite BAMs
    // First, count expected BAMs per sample for non-blocking groupTuple op
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

    // Now, group with expected size then sort into tumor and normal channels
    // channel: [ meta_group, [aln, ...], [idx, ...] ]
    ch_alns_united = ch_sample_fastq_counts
        // channel: [ [ meta_group, count ], [ meta_group, aln, idx ] ]
        .cross(
            // First element to match meta_group above for `cross`
            channel.topic('bwamem2_align_bam').map { meta_bwamem2, aln, idx -> [[key: meta_bwamem2.key, sample_type: meta_bwamem2.sample_type], aln, idx] }
        )
        .map { count_tuple, inputs_tuple ->
            def group_size = count_tuple[1]
            def (meta_group, aln, idx) = inputs_tuple

            return tuple(groupKey(meta_group, group_size), aln, idx)
        }
        .groupTuple()
        .branch { meta_group, alns, idxs ->
            assert ['tumor', 'normal', 'donor'].contains(meta_group.sample_type)
            tumor: meta_group.sample_type == 'tumor'
            normal: meta_group.sample_type == 'normal'
            donor: meta_group.sample_type == 'donor'
            placeholder: true
        }

    //
    // STEP: Handle outputs
    //
    // Set outputs, restoring original meta
    // channel: [ meta, [aln, ...], [idx, ...] ]
    ch_outputs_tumor = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_alns_united.tumor, ch_inputs),
            ch_inputs_tumor_sorted.skip.unique().map { meta -> [meta, [], []] },
        )

    // channel: [ meta, [aln, ...], [idx, ...] ]
    ch_outputs_normal = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_alns_united.normal, ch_inputs),
            ch_inputs_normal_sorted.skip.unique().map { meta -> [meta, [], []] },
        )

    // channel: [ meta, [aln, ...], [idx, ...] ]
    ch_outputs_donor = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_alns_united.donor, ch_inputs),
            ch_inputs_donor_sorted.skip.unique().map { meta -> [meta, [], []] },
        )

    emit:
    tumor  = ch_outputs_tumor  // channel: [ meta, [aln, ...], [idx, ...] ]
    normal = ch_outputs_normal // channel: [ meta, [aln, ...], [idx, ...] ]
    donor  = ch_outputs_donor  // channel: [ meta, [aln, ...], [idx, ...] ]
}
