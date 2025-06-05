//
// Align RNA reads
//

import Constants
import Utils

include { GATK4_MARKDUPLICATES } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { SAMBAMBA_MERGE       } from '../../../modules/local/sambamba/merge/main'
include { SAMTOOLS_SORT        } from '../../../modules/nf-core/samtools/sort/main'
include { STAR_ALIGN           } from '../../../modules/local/star/align/main'

workflow READ_ALIGNMENT_RNA {
    take:
    // Sample data
    ch_inputs         // channel: [mandatory] [ meta ]

    // Reference data
    genome_star_index // channel: [mandatory] /path/to/genome_star_index/

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Sort inputs
    // channel: [ meta ]
    ch_inputs_sorted = ch_inputs
        .branch { meta ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAM_RNA_TUMOR)
            runnable: Utils.hasTumorRnaFastq(meta) && !has_existing
            skip: true
        }

    // Create FASTQ input channel
    // channel: [ meta_fastq, fastq_fwd, fastq_rev ]
    ch_fastq_inputs = ch_inputs_sorted.runnable
        .flatMap { meta ->
            def meta_sample = Utils.getTumorRnaSample(meta)
            meta_sample
                .getAt(Constants.FileType.FASTQ)
                .collect { key, fps ->
                    def (library_id, lane) = key

                    def meta_fastq = [
                        key: meta.group_id,
                        id: "${meta.group_id}_${meta_sample.sample_id}",
                        sample_id: meta_sample.sample_id,
                        library_id: library_id,
                        lane: lane,
                    ]

                    return [meta_fastq, fps['fwd'], fps['rev']]
                }
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

    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

    //
    // MODULE: SAMtools sort
    //
    // Create process input channel
    // channel: [ meta_sort, bam ]
    ch_sort_inputs = STAR_ALIGN.out.bam
        .map { meta_star, bam ->
            def meta_sort = [
                *:meta_star,
                prefix: meta_star.read_group,
            ]

            return [meta_sort, bam]
        }

    // Run process
    SAMTOOLS_SORT(
        ch_sort_inputs,
    )

    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

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
    // channel: [ meta_group, [bam, ...] ]
    ch_bams_united = ch_sample_fastq_counts
        .cross(
            // First element to match meta_count above for `cross`
            SAMTOOLS_SORT.out.bam.map { meta_star, bam -> [[key: meta_star.key], bam] }
        )
        .map { count_tuple, bam_tuple ->

            def group_size = count_tuple[1]
            def (meta_bam, bam) = bam_tuple

            def meta_group = [
                *:meta_bam,
            ]

            return tuple(groupKey(meta_group, group_size), bam)
        }
        .groupTuple()

    // Sort into merge-eligible BAMs (at least two BAMs required)
    // channel: runnable: [ meta_group, [bam, ...] ]
    // channel: skip: [ meta_group, bam ]
    ch_bams_united_sorted = ch_bams_united
        .branch { meta_group, bams ->
            runnable: bams.size() > 1
            skip: true
                return [meta_group, bams[0]]
        }

    // Create process input channel
    // channel: [ meta_merge, [bams, ...] ]
    ch_merge_inputs = WorkflowOncoanalyser.restoreMeta(ch_bams_united_sorted.runnable, ch_inputs)
        .map { meta, bams ->
            def meta_merge = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorRnaSampleName(meta),
            ]
            return [meta_merge, bams]
        }

    // Run process
    SAMBAMBA_MERGE(
        ch_merge_inputs,
    )

    ch_versions = ch_versions.mix(SAMBAMBA_MERGE.out.versions)

    //
    // MODULE: GATK4 markduplicates
    //
    // Create process input channel
    // channel: [ meta_markdups, bam ]
    ch_markdups_inputs = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(SAMBAMBA_MERGE.out.bam, ch_inputs),
            WorkflowOncoanalyser.restoreMeta(ch_bams_united_sorted.skip, ch_inputs),
        )
        .map { meta, bam ->
            def meta_markdups = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorRnaSampleName(meta),
            ]
            return [meta_markdups, bam]
        }

    // Run process
    GATK4_MARKDUPLICATES(
        ch_markdups_inputs,
        [],
        [],
    )

    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions)

    // Combine BAMs and BAIs
    // channel: [ meta, bam, bai ]
    ch_bams_ready = WorkflowOncoanalyser.groupByMeta(
        WorkflowOncoanalyser.restoreMeta(GATK4_MARKDUPLICATES.out.bam, ch_inputs),
        WorkflowOncoanalyser.restoreMeta(GATK4_MARKDUPLICATES.out.bai, ch_inputs),
    )

    // Set outputs
    // channel: [ meta, bam, bai ]
    ch_bam_out = Channel.empty()
        .mix(
            ch_bams_ready,
            ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
        )

    emit:
    rna_tumor = ch_bam_out  // channel: [ meta, bam, bai ]

    versions  = ch_versions // channel: [ versions.yml ]
}
