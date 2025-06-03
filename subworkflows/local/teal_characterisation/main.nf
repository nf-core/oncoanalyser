//
// TEAL performs characterisation of telomeric features and rearrangements
//

import Constants
import Utils

include { TEAL } from '../../../modules/local/teal/main'

workflow TEAL_CHARACTERISATION {
    take:
    // Sample data
    ch_inputs         // channel: [mandatory] [ meta ]
    ch_tumor_bam      // channel: [mandatory] [ meta, bam, bai ]
    ch_normal_bam     // channel: [mandatory] [ meta, bam, bai ]
    ch_tumor_metrics  // channel: [mandatory] [ meta, metrics_dir ]
    ch_normal_metrics // channel: [mandatory] [ meta, metrics_dir ]
    ch_cobalt         // channel: [mandatory] [ meta, cobalt_dir ]
    ch_purple         // channel: [mandatory] [ meta, purple_dir ]

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources and sort
    // channel: runnable: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai, tumor_metrics_dir, normal_metrics_dir, cobalt_dir, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_tumor_bam,
        ch_normal_bam,
        ch_tumor_metrics,
        ch_normal_metrics,
        ch_cobalt,
        ch_purple,
    )
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, tumor_metrics_dir, normal_metrics_dir, cobalt_dir, purple_dir ->
            return [
                meta,
                Utils.selectCurrentOrExisting(tumor_bam, meta, Constants.INPUT.BAM_REDUX_DNA_TUMOR),
                tumor_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_TUMOR),
                Utils.selectCurrentOrExisting(normal_bam, meta, Constants.INPUT.BAM_REDUX_DNA_NORMAL),
                normal_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_NORMAL),
                Utils.selectCurrentOrExisting(tumor_metrics_dir, meta, Constants.INPUT.BAMTOOLS_DIR_TUMOR),
                Utils.selectCurrentOrExisting(normal_metrics_dir, meta, Constants.INPUT.BAMTOOLS_DIR_NORMAL),
                Utils.selectCurrentOrExisting(cobalt_dir, meta, Constants.INPUT.COBALT_DIR),
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
            ]
        }
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, tumor_metrics_dir, normal_metrics_dir, cobalt_dir, purple_dir ->

            def has_tumor = tumor_bam && tumor_metrics_dir && purple_dir
            def has_normal = normal_bam && normal_metrics_dir
            def has_input = tumor_bam || normal_bam

            runnable: has_input && cobalt_dir
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_teal, tumor_bam, tumor_bai, normal_bam, normal_bai, tumor_metrics_dir, normal_metrics_dir, cobalt_dir, purple_dir ]
    ch_teal_inputs = ch_inputs_sorted.runnable
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, tumor_metrics_dir, normal_metrics_dir, cobalt_dir, purple_dir ->

            def meta_teal = [
                key: meta.group_id,
                id: meta.group_id,
            ]

            if (tumor_bam) {
                meta_teal.tumor_id = Utils.getTumorDnaSampleName(meta)
            }

            if (normal_bam) {
                meta_teal.normal_id = Utils.getNormalDnaSampleName(meta)
            }

            return [meta_teal, tumor_bam, tumor_bai, normal_bam, normal_bai, tumor_metrics_dir, normal_metrics_dir, cobalt_dir, purple_dir]
        }

    // Run process
    TEAL(
        ch_teal_inputs,
    )

    ch_versions = ch_versions.mix(TEAL.out.versions)

    emit:
    versions  = ch_versions // channel: [ versions.yml ]
}
