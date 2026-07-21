//
// COBALT calculates read ratios between tumor and normal samples
//

import Constants
import Utils

include { COBALT } from '../../../modules/local/cobalt/run/main'

workflow COBALT_PROFILING {
    take:
    // Sample data
    ch_inputs                    // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor           // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal          // channel: [mandatory] [ meta, redux_dir ]

    // Reference data
    genome_version               // channel: [mandatory] genome version
    gc_profile                   // channel: [mandatory] /path/to/gc_profile
    diploid_bed                  // channel: [optional]  /path/to/diploid_bed
    target_regions_normalisation // channel: [optional]  /path/to/target_regions_normalisation
    targeted_mode                // boolean: [mandatory] Set targeted mode
    purity_estimate_mode         // boolean: [mandatory] Set purity estimate mode

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources then sort
    // NOTE(SW): germline mode is not currently supported
    // channel: runnable: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_redux_dir_tumor,
        ch_redux_dir_normal,
    )
        .map { meta, redux_dir_tumor, redux_dir_normal ->

            def redux_dir_tumor_selected = Utils.selectCurrentOrExisting(redux_dir_tumor, meta, Constants.INPUT.REDUX_DIR_TUMOR)
            def redux_dir_normal_selected = Utils.selectCurrentOrExisting(redux_dir_normal, meta, Constants.INPUT.REDUX_DIR_NORMAL)

            def (tumor_bam, tumor_bai) = Utils.getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)
            def (normal_bam, normal_bai) = Utils.getNormalReduxDirAlignment(meta, redux_dir_normal_selected)

            return [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]

        }
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.COBALT_DIR)
            runnable_tn: tumor_bam && normal_bam && ! has_existing
            runnable_to: tumor_bam && ! has_existing
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
    // channel: sample_data: [ meta_cobalt, tumor_bam, tumor_bai, normal_bam, normal_bai ]
    // channel: diploid_bed: [ diploid_bed ]
    ch_cobalt_inputs = ch_inputs_runnable
        .multiMap { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, diploid_bed ->

            def tumor_id = purity_estimate_mode
            if (purity_estimate_mode) {
                tumor_id = Utils.getTumorDnaSampleName(meta, primary: false)
            } else {
                tumor_id = Utils.getTumorDnaSampleName(meta, primary: true)
            }

            def meta_cobalt = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: tumor_id,
            ]

            if (normal_bam) {
                meta_cobalt.normal_id = Utils.getNormalDnaSampleName(meta)
            }

            sample_data: [meta_cobalt, tumor_bam, tumor_bai, normal_bam, normal_bai]
            diploid_bed: diploid_bed
        }

    // Run process
    COBALT(
        ch_cobalt_inputs.sample_data,
        genome_version,
        gc_profile,
        ch_cobalt_inputs.diploid_bed,
        target_regions_normalisation,
        targeted_mode,
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
