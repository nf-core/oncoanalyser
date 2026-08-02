//
// Align RNA reads
//

include { GATK4_MARKDUPLICATES } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { SAMBAMBA_MERGE       } from '../../../modules/local/sambamba/merge/main'
include { SAMTOOLS_SORT        } from '../../../modules/nf-core/samtools/sort/main'
include { STAR_ALIGN           } from '../../../modules/local/star/align/main'

workflow READ_ALIGNMENT_RNA {
    take:
    // Sample data
    ch_inputs         // channel: [mandatory] [ meta ]
    ch_fastq          // channel: [mandatory] [ meta, fastq_info, fastq_fwd, fastq_rev ]

    // Reference data
    genome_star_index // channel: [mandatory] /path/to/genome_star_index/

    main:
    //
    // STEP: Handle inputs
    //
    // Sort inputs
    // runnable: channel: [ meta, fastq_info, fastq_fwd, fastq_rev ]
    // skip: channel: [ meta ]
    ch_inputs_sorted = ch_fastq
        .branch { meta, fastq_info, fastq_fwd, fastq_rev ->
            def has_inputs = fastq_fwd && fastq_rev
            runnable: has_inputs
            skip: true
              return meta
        }

    // Create FASTQ input channel
    // channel: [ meta_fastq, fastq_fwd, fastq_rev ]
    ch_fastq_inputs = ch_inputs_sorted.runnable
        .map { meta, fastq_info, fastq_fwd, fastq_rev ->

            def meta_fastq = [
                key: meta.group_id,
                id: "${meta.group_id}_${fastq_info.sample_id}",
                sample_id: fastq_info.sample_id,
                library_id: fastq_info.library_id,
                lane: fastq_info.lane,
            ]

            return [meta_fastq, fastq_fwd, fastq_rev]
        }

    //
    // MODULE: STAR alignment
    //
    // Create process input channel
    // channel: [ meta_star, fastq_fwd, fastq_rev ]
    ch_star_inputs = ch_fastq_inputs
        .map { meta_fastq, fastq_fwd, fastq_rev ->
            def meta_star = [
                *:meta_fastq,
                read_group: "${meta_fastq.sample_id}.${meta_fastq.library_id}.${meta_fastq.lane}",
            ]

            return [meta_star, fastq_fwd, fastq_rev]
        }

    // Run process
    STAR_ALIGN(
        ch_star_inputs,
        genome_star_index,
    )

    //
    // MODULE: SAMtools sort
    //
    // Create process input channel
    // channel: [ meta_sort, aln ]
    ch_sort_inputs = channel.topic('star_align_bam')
        .map { meta_star, aln ->
            def meta_sort = meta_star + [prefix: meta_star.sample_id]
            return [meta_sort, aln]
        }

    // Run process
    SAMTOOLS_SORT(
        ch_sort_inputs,
    )

    //
    // MODULE: Sambamba merge
    //
    // Reunite BAMs
    // First, count expected BAMs per sample for non-blocking groupTuple op
    // channel: [ meta_count, group_size ]
    ch_sample_fastq_counts = ch_star_inputs
        .map { meta_star, reads_fwd, reads_rev ->
            def meta_count = [key: meta_star.key]
            return [meta_count, meta_star]
        }
        .groupTuple()
        .map { meta_count, meta_stars -> return [meta_count, meta_stars.size()] }

    // Now, group with expected size then sort into tumor and normal channels
    // channel: [ meta_group, [aln, ...] ]
    ch_alns_united = ch_sample_fastq_counts
        .cross(
            // First element to match meta_count above for `cross`
            SAMTOOLS_SORT.out.bam.map { meta_star, aln -> [[key: meta_star.key], aln] }
        )
        .map { count_tuple, aln_tuple ->

            def group_size = count_tuple[1]
            def (meta_aln, aln) = aln_tuple

            def meta_group = [
                *:meta_aln,
            ]

            return tuple(groupKey(meta_group, group_size), aln)
        }
        .groupTuple()

    // Sort into merge-eligible BAMs (at least two BAMs required)
    // channel: runnable: [ meta_group, [aln, ...] ]
    // channel: skip: [ meta_group, aln ]
    ch_alns_united_sorted = ch_alns_united
        .branch { meta_group, alns ->
            runnable: alns.size() > 1
            skip: true
                return [meta_group, alns[0]]
        }

    // Create process input channel
    // channel: [ meta_merge, [alns, ...] ]
    ch_merge_inputs = WorkflowOncoanalyser.restoreMeta(ch_alns_united_sorted.runnable, ch_inputs)
        .map { meta, alns ->
            def meta_merge = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorRnaSampleName(meta),
            ]
            return [meta_merge, alns]
        }

    // Run process
    SAMBAMBA_MERGE(
        ch_merge_inputs,
    )

    //
    // MODULE: GATK4 markduplicates
    //
    // Create process input channel
    // channel: [ meta_markdups, aln ]
    ch_markdups_inputs = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('sambamba_merge_bam'), ch_inputs),
            WorkflowOncoanalyser.restoreMeta(ch_alns_united_sorted.skip, ch_inputs),
        )
        .map { meta, aln ->
            def meta_markdups = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorRnaSampleName(meta),
            ]
            return [meta_markdups, aln]
        }

    // Run process
    GATK4_MARKDUPLICATES(
        ch_markdups_inputs,
        [],
        [],
    )

    //
    // STEP: Handle outputs
    //
    // Combine BAMs and BAIs
    // channel: [ meta, aln, idx ]
    ch_alns_ready = WorkflowOncoanalyser.groupByMeta(
        WorkflowOncoanalyser.restoreMeta(channel.topic('gatk4_markduplicates_bam'), ch_inputs),
        WorkflowOncoanalyser.restoreMeta(channel.topic('gatk4_markduplicates_bai'), ch_inputs),
    )

    // Combine STAR log with QC and MarkDuplicates metrics
    // channel: [ meta, star_log, md_metrics ]
    ch_qc_files_ready = WorkflowOncoanalyser.groupByMeta(
        WorkflowOncoanalyser.restoreMeta(channel.topic('star_align_qc_log'), ch_inputs),
        WorkflowOncoanalyser.restoreMeta(channel.topic('gatk4_markduplicates_metrics'), ch_inputs),
    )

    // Set outputs
    // channel: [ meta, aln, idx ]
    ch_outputs_aln = channel.empty()
        .mix(
            ch_alns_ready,
            ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
        )

    // channel: [ meta, star_log, md_metrics ]
    ch_outputs_qc_files = channel.empty()
        .mix(
            ch_qc_files_ready,
            ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
        )

    emit:
    tumor    = ch_outputs_aln      // channel: [ meta, aln, idx ]
    qc_files = ch_outputs_qc_files // channel: [ meta, star_log, md_metrics ]
}
