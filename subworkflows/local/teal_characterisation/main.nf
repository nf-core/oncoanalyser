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
    ch_inputs              // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor     // channel: [mandatory] [ meta, bam, bai ]
    ch_redux_dir_normal    // channel: [mandatory] [ meta, bam, bai ]
    ch_bamtools_dir_tumor  // channel: [mandatory] [ meta, bamtools_dir ]
    ch_bamtools_dir_normal // channel: [mandatory] [ meta, bamtools_dir ]
    ch_cobalt_dir          // channel: [mandatory] [ meta, cobalt_dir ]
    ch_purple_dir          // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_version         // channel: [mandatory] genome version

    // Params
    sequencing_platform    // string:  [mandatory] sequencing platform

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    //
    // STEP: Handle inputs
    //
    // Select input sources then sort
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
            runnable: tumor_bam || normal_bam
            skip: true
                return meta
        }

    //
    // MODULE: TEAL prep
    //
    // Create process input channel
    // channel: [ meta_teal, tumor_bam, tumor_bai, normal_bam, normal_bai ]
    ch_teal_prep_inputs = ch_inputs_sorted.runnable
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->

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

            return [meta_teal, tumor_bam, tumor_bai, normal_bam, normal_bai]
        }

    // Run process
    TEAL_PREP(
        ch_teal_prep_inputs,
        genome_version,
        sequencing_platform,
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
    // MODULE: TEAL pipeline
    //
    // Select input sources then sort
    // channel: runnable: [ meta, teal_bam_tumor, teal_bai_tumor, teal_bam_normal, teal_bai_normal, bamtools_dir_tumor, bamtools_dir_normal, cobalt_dir, purple_dir ]
    // channel: skip: [ meta ]
    ch_teal_pipeline_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_tumor_teal_bam,
        ch_normal_teal_bam,
        ch_bamtools_dir_tumor,
        ch_bamtools_dir_normal,
        ch_cobalt_dir,
        ch_purple_dir,
    )
        .map { meta, teal_bam_tumor, teal_bai_tumor, teal_bam_normal, teal_bai_normal, bamtools_dir_tumor, bamtools_dir_normal, cobalt_dir, purple_dir ->
            return [
                meta,
                teal_bam_tumor,
                teal_bai_tumor,
                teal_bam_normal,
                teal_bai_normal,
                Utils.selectCurrentOrExisting(bamtools_dir_tumor, meta, Constants.INPUT.BAMTOOLS_DIR_TUMOR),
                Utils.selectCurrentOrExisting(bamtools_dir_normal, meta, Constants.INPUT.BAMTOOLS_DIR_NORMAL),
                Utils.selectCurrentOrExisting(cobalt_dir, meta, Constants.INPUT.COBALT_DIR),
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
            ]
        }
        .branch { meta, teal_bam_tumor, teal_bai_tumor, teal_bam_normal, teal_bai_normal, bamtools_dir_tumor, bamtools_dir_normal, cobalt_dir, purple_dir ->

            def has_tumor = teal_bam_tumor && bamtools_dir_tumor && purple_dir
            def has_normal = teal_bam_normal && bamtools_dir_normal

            runnable: (has_tumor || has_normal) && cobalt_dir
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_teal, teal_bam_tumor, teal_bai_tumor, teal_bam_normal, teal_bai_normal, bamtools_dir_tumor, bamtools_dir_normal, cobalt_dir, purple_dir ]
    ch_teal_pipeline_inputs = ch_teal_pipeline_inputs_sorted.runnable
        .map { meta, teal_bam_tumor, teal_bai_tumor, teal_bam_normal, teal_bai_normal, bamtools_dir_tumor, bamtools_dir_normal, cobalt_dir, purple_dir ->

            def meta_teal = [
                key: meta.group_id,
                id: meta.group_id,
            ]

            if (teal_bam_tumor) {
                meta_teal.tumor_id = Utils.getTumorDnaSampleName(meta)
            }

            if (teal_bam_normal) {
                meta_teal.normal_id = Utils.getNormalDnaSampleName(meta)
            }

            return [meta_teal, teal_bam_tumor, teal_bai_tumor, teal_bam_normal, teal_bai_normal, bamtools_dir_tumor, bamtools_dir_normal, cobalt_dir, purple_dir]
        }

    // Run process
    TEAL_PIPELINE(
        ch_teal_pipeline_inputs,
        genome_version,
        sequencing_platform,
    )

    ch_versions = ch_versions.mix(TEAL_PIPELINE.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
