//
// Align RNA reads from an existing BAM by streaming through samtools fastq into STAR
//

import Constants
import Utils

include { GATK4_MARKDUPLICATES  } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { SAMTOOLS_SORT         } from '../../../modules/nf-core/samtools/sort/main'
include { STAR_ALIGN_FROM_BAM   } from '../../../modules/local/star/align_from_bam/main'

workflow READ_ALIGNMENT_RNA_FROM_BAM {
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
            runnable: Utils.hasTumorRnaBam(meta)
            skip: true
        }

    // Create BAM/CRAM input channel for STAR_ALIGN_FROM_BAM
    // channel: [ meta_star, bam ]
    ch_star_inputs = ch_inputs_sorted.runnable
        .map { meta ->
            def meta_sample = Utils.getTumorRnaSample(meta)
            def meta_star = [
                key:        meta.group_id,
                id:         "${meta.group_id}_${meta_sample.sample_id}",
                sample_id:  meta_sample.sample_id,
                read_group: "${meta_sample.sample_id}.realign",
            ]
            return [meta_star, Utils.getTumorRnaBam(meta)]
        }

    STAR_ALIGN_FROM_BAM(
        ch_star_inputs,
        genome_star_index,
    )

    ch_versions = ch_versions.mix(STAR_ALIGN_FROM_BAM.out.versions)

    // Sort the unsorted BAM output from STAR
    ch_sort_inputs = STAR_ALIGN_FROM_BAM.out.bam
        .map { meta_star, bam ->
            [[ *:meta_star, prefix: meta_star.read_group ], bam]
        }

    SAMTOOLS_SORT(ch_sort_inputs)

    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    // Markduplicates (single BAM per sample, no merge needed)
    ch_markdups_inputs = WorkflowOncoanalyser.restoreMeta(SAMTOOLS_SORT.out.bam, ch_inputs)
        .map { meta, bam ->
            def meta_markdups = [
                key:       meta.group_id,
                id:        meta.group_id,
                sample_id: Utils.getTumorRnaSampleName(meta),
            ]
            return [meta_markdups, bam]
        }

    GATK4_MARKDUPLICATES(ch_markdups_inputs, [], [])

    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions)

    ch_bams_ready = WorkflowOncoanalyser.groupByMeta(
        WorkflowOncoanalyser.restoreMeta(GATK4_MARKDUPLICATES.out.bam, ch_inputs),
        WorkflowOncoanalyser.restoreMeta(GATK4_MARKDUPLICATES.out.bai, ch_inputs),
    )

    ch_bam_out = Channel.empty()
        .mix(
            ch_bams_ready,
            ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
        )

    emit:
    rna_tumor = ch_bam_out   // channel: [ meta, bam, bai ]

    versions  = ch_versions  // channel: [ versions.yml ]
}
