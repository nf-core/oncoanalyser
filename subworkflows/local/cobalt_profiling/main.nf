//
// COBALT calculates read ratios between tumor and normal samples
//

import Constants
import Utils

include { COBALT } from '../../../modules/local/cobalt/main'

workflow COBALT_PROFILING {
    take:
    // Sample data
    ch_inputs                   // channel: [mandatory] [ meta ]
    ch_tumor_bam                // channel: [mandatory] [ meta, bam, bai ]
    ch_normal_bam               // channel: [mandatory] [ meta, bam, bai ]

    // Reference data
    gc_profile                  // channel: [mandatory] /path/to/gc_profile
    diploid_bed                 // channel: [optional]  /path/to/diploid_bed
    target_region_normalisation // channel: [optional]  /path/to/target_region_normalisation

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources and sort
    // NOTE(SW): germline mode is not currently supported
    // channel: runnable: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_tumor_bam,
        ch_normal_bam,
    )
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            return [
                meta,
                Utils.selectCurrentOrExisting(tumor_bam, meta, Constants.INPUT.BAM_REDUX_DNA_TUMOR),
                tumor_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_TUMOR),
                Utils.selectCurrentOrExisting(normal_bam, meta, Constants.INPUT.BAM_REDUX_DNA_NORMAL),
                normal_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_NORMAL),
            ]
        }
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.COBALT_DIR)
            runnable_tn: tumor_bam && normal_bam && !has_existing
            runnable_to: tumor_bam && !has_existing
            skip: true
                return meta
        }

    // First set diploid BED input for tumor/normal and tumor only samples
    // NOTE(SW): since the diploid BED is provided as a channel, I seem to be only able to include via channel ops
    // channel: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai, diploid_bed ]
    ch_inputs_runnable = Channel.empty()
        .mix(
            ch_inputs_sorted.runnable_tn.map { [*it, []] },
            ch_inputs_sorted.runnable_to.combine(diploid_bed),
        )

    // Create process input channel
    // channel: sample_data: [ meta_cobalt, tumor_bam, normal_bam, tumor_bai, normal_bai ]
    // channel: diploid_bed: [ diploid_bed ]
    ch_cobalt_inputs = ch_inputs_runnable
        .multiMap { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, diploid_bed ->

            def meta_cobalt = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Utils.getTumorDnaSampleName(meta),
            ]

            if (normal_bam) {
                meta_cobalt.normal_id = Utils.getNormalDnaSampleName(meta)
            }

            sample_data: [meta_cobalt, tumor_bam, normal_bam, tumor_bai, normal_bai]
            diploid_bed: diploid_bed
        }

    // Run process
    COBALT(
        ch_cobalt_inputs.sample_data,
        gc_profile,
        ch_cobalt_inputs.diploid_bed,
        target_region_normalisation,
    )

    ch_versions = ch_versions.mix(COBALT.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, cobalt_dir ]
    ch_outputs = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(COBALT.out.cobalt_dir, ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    cobalt_dir = ch_outputs  // channel: [ meta, cobalt_dir ]

    versions   = ch_versions // channel: [ versions.yml ]
}
