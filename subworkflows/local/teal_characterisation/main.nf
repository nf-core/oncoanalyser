//
// TEAL performs characterisation of telomeric features and rearrangements
//

import Constants
import Utils

include { TEAL_PREP     } from '../../../modules/local/teal/prep/main'
include { TEAL_PIPELINE } from '../../../modules/local/teal/pipeline/main'

workflow TEAL_CHARACTERISATION {
    take:
    // Sample data
    ch_inputs         // channel: [mandatory] [ meta ]
    ch_tumor_bam      // channel: [mandatory] [ meta, bam, bai ]
    ch_normal_bam     // channel: [mandatory] [ meta, bam, bai ]
    ch_tumor_metrics  // channel: [mandatory] [ meta, metrics_dir ]
    ch_normal_metrics // channel: [mandatory] [ meta, metrics_dir ]
    ch_cobalt_dir     // channel: [mandatory] [ meta, cobalt_dir ]
    ch_purple_dir     // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_version    // channel: [mandatory] genome version

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    //
    // MODULE: Teal prep
    //
    // Make preliminary BAM to be used as input for the main Teal pipeline. This prep step takes 90% of Teal's runtime,
    // so it is split from the main pipeline so that we don't have to wait for Purple to finish before we start
    // running Teal
    //

    // Select input sources and sort
    // channel: runnable: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai ]
    // channel: skip: [ meta ]
    ch_teal_prep_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
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
            runnable: tumor_bam || normal_bam
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_teal, tumor_bam, tumor_bai, normal_bam, normal_bai ]
    ch_teal_prep_inputs = ch_teal_prep_inputs_sorted.runnable
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->

            def meta_teal = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: tumor_bam ? Utils.getTumorDnaSampleName(meta) : null,
                normal_id: normal_bam ? Utils.getNormalDnaSampleName(meta) : null,
            ]

            return [meta_teal, tumor_bam, tumor_bai, normal_bam, normal_bai]
        }

    // Run process
    TEAL_PREP(
        ch_teal_prep_inputs,
        genome_version,
    )

    ch_versions = ch_versions.mix(TEAL_PREP.out.versions)

    // Flatten TEAL_PREP output
    // channel: [ meta, teal_bam, teal_bai ]
    ch_tumor_teal_bam = WorkflowOncoanalyser.restoreMeta(TEAL_PREP.out.tumor_bam, ch_inputs)
        .map { meta, bam_bai -> [meta, *bam_bai] }

    ch_normal_teal_bam_placeholder = WorkflowOncoanalyser.restoreMeta(
        ch_teal_prep_inputs
            .filter { it[0].normal_id == null } // Only populate placeholder channel if normal sample is missing
            .map { [ it[0], [], [] ] },
        ch_inputs
    )

    ch_normal_teal_bam = WorkflowOncoanalyser.restoreMeta(TEAL_PREP.out.normal_bam, ch_inputs)
        .map { meta, bam_bai -> [meta, *bam_bai] }
        .mix(ch_normal_teal_bam_placeholder)

    //
    // MODULE: Teal pipeline
    //

    // Select input sources and sort
    // channel: runnable: [ meta, tumor_teal_bam, tumor_teal_bai, normal_teal_bam, normal_teal_bai, tumor_metrics_dir, normal_metrics_dir, cobalt_dir, purple_dir ]
    // channel: skip: [ meta ]
    ch_teal_pipeline_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_tumor_teal_bam,
        ch_normal_teal_bam,
        ch_tumor_metrics,
        ch_normal_metrics,
        ch_cobalt_dir,
        ch_purple_dir,
    )
        .map { meta, tumor_teal_bam, tumor_teal_bai, normal_teal_bam, normal_teal_bai, tumor_metrics_dir, normal_metrics_dir, cobalt_dir, purple_dir ->

            return [
                meta,
                tumor_teal_bam, tumor_teal_bai,
                normal_teal_bam, normal_teal_bai,
                Utils.selectCurrentOrExisting(tumor_metrics_dir, meta, Constants.INPUT.BAMTOOLS_DIR_TUMOR),
                Utils.selectCurrentOrExisting(normal_metrics_dir, meta, Constants.INPUT.BAMTOOLS_DIR_NORMAL),
                Utils.selectCurrentOrExisting(cobalt_dir, meta, Constants.INPUT.COBALT_DIR),
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
            ]
        }
        .branch { meta, tumor_teal_bam, tumor_teal_bai, normal_teal_bam, normal_teal_bai, tumor_metrics_dir, normal_metrics_dir, cobalt_dir, purple_dir ->

            def has_tumor = tumor_teal_bam && tumor_metrics_dir && purple_dir
            def has_normal = normal_teal_bam && normal_metrics_dir

            runnable: (has_tumor || has_normal) && cobalt_dir
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_teal, tumor_teal_bam, tumor_teal_bai, normal_teal_bam, normal_teal_bai, tumor_metrics_dir, normal_metrics_dir, cobalt_dir, purple_dir ]
    ch_teal_pipeline_inputs = ch_teal_pipeline_inputs_sorted.runnable
        .map { meta, tumor_teal_bam, tumor_teal_bai, normal_teal_bam, normal_teal_bai, tumor_metrics_dir, normal_metrics_dir, cobalt_dir, purple_dir ->

            def meta_teal = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: tumor_teal_bam ? Utils.getTumorDnaSampleName(meta) : null,
                normal_id: normal_teal_bam ? Utils.getNormalDnaSampleName(meta) : null,
            ]

            return [ meta_teal, tumor_teal_bam, tumor_teal_bai, normal_teal_bam, normal_teal_bai, tumor_metrics_dir, normal_metrics_dir, cobalt_dir, purple_dir ]
        }


    TEAL_PIPELINE(
        ch_teal_pipeline_inputs,
        genome_version,
    )

    ch_versions = ch_versions.mix(TEAL_PIPELINE.out.versions)


    emit:
    versions  = ch_versions // channel: [ versions.yml ]
}
