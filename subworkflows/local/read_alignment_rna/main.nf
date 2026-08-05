//
// Align RNA reads
//

include { GATK4_MARKDUPLICATES } from '../../../modules/nf-core/gatk4/markduplicates/main'
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

            // NOTE(SW): initial map sets defaults and conventional ordering of selected fields, merging then overwrites / adds while preserving order
            def rg_id = [fastq_info.sample_id, fastq_info.library_id, fastq_info.lane, fastq_info.flowcell].findAll().join('.')
            def rg_entries = [ID: rg_id, SM: fastq_info.sample_id, LB: fastq_info.library_id] + fastq_info.rg_fields
            def rg_line = rg_entries.collect { k, v -> "'${k}:${v}'" }.join(' ')

            def meta_fastq = [
                key: meta.group_id,
                id: "${meta.group_id}_${fastq_info.sample_id}",
                sample_id: fastq_info.sample_id,
                rg_line: rg_line,
            ]

            return [meta_fastq, fastq_fwd, fastq_rev]
        }

    //
    // MODULE: STAR alignment
    //
    // Create process input channel
    // First, count expected FASTQ pairs per sample for non-blocking groupTuple op
    // channel: [ meta_group, group_size ]
    ch_fastq_counts = ch_fastq_inputs
        .map { meta_fastq, _fastq_fwd, _fastq_rev ->
            def meta_group = [key: meta_fastq.key]
            return [meta_group, meta_fastq]
        }
        .groupTuple()
        .map { meta_group, meta_fastqs -> return [meta_group, meta_fastqs.size()] }

    // Now, group with expected size to proceed without blocking
    // channel: [ meta_star, [ rg_line, ... ], [ fastq_fwd, ... ], [ fastq_rev, ... ] ]
    ch_star_inputs = ch_fastq_counts
        // channel: [ [ meta_group, count ], [ meta_group, fastq_fwd, fastq_rev ] ]
        .cross(
            // First element to match meta_group above for `cross`
            ch_fastq_inputs.map { meta_fastq, fastq_fwd, fastq_rev -> [[key: meta_fastq.key], meta_fastq, fastq_fwd, fastq_rev] }
        )
        // channel: [ GroupKey(meta_group, size), rg_line, fastq_fwd, fastq_rev ]
        .map { count_tuple, inputs_tuple ->
            def group_size = count_tuple[1]
            def (meta_group, meta_fastq, fastq_fwd, fastq_rev) = inputs_tuple

            def meta_star = [
                key: meta_fastq.key,
                id: meta_fastq.id,
                sample_id: meta_fastq.sample_id,
            ]

            return tuple(groupKey(meta_star, group_size), meta_fastq.rg_line, fastq_fwd, fastq_rev)
        }
        // channel: [ meta_star, rg_lines, fastq_fwds, fastq_revs ]
        .groupTuple()
        // Sort for stablised inputs
        .map { meta_star, rg_lines, fastq_fwds, fastq_revs ->
            def entries_sorted = [rg_lines, fastq_fwds, fastq_revs]
                // [ [ rg_line, fastq_fwd, fastq_rev ], ... ]
                .transpose()
                .sort()
                // [ [ rg_line, ... ], [ fastq_fwd, ... ], [ fastq_rev, ... ] ]
                .transpose()

            // Unpack cleanly
            def (rg_lines_sorted, fastq_fwds_sorted, fastq_revs_sorted) = entries_sorted

            return [meta_star, rg_lines_sorted, fastq_fwds_sorted, fastq_revs_sorted]
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
    // MODULE: GATK4 markduplicates
    //
    // Create process input channel
    // channel: [ meta_markdups, aln ]
    ch_markdups_inputs = WorkflowOncoanalyser.restoreMeta(channel.topic('samtools_sort_bam'), ch_inputs)
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
