//
// LINX plotting visualises clusters structural variants
//

import Constants
import Utils

include { LINXREPORT      } from '../../../modules/local/linxreport/main'
include { LINX_VISUALISER } from '../../../modules/local/linx/visualiser/main'

workflow LINX_PLOTTING {
    take:
    // Sample data
    ch_inputs              // channel: [mandatory] [ meta ]
    ch_annotations         // channel: [mandatory] [ meta, annotation_dir ]

    // Reference data
    genome_version         // channel: [mandatory] genome version
    ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/

    main:
    // Channel for versions.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources and sort
    // channel: runnable: [ meta, annotation_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_annotations
        .map { meta, annotation_dir ->
            return [
                meta,
                Utils.selectCurrentOrExisting(annotation_dir, meta, Constants.INPUT.LINX_ANNO_DIR_TUMOR),
            ]
        }
        .branch { meta, annotation_dir ->

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.LINX_PLOT_DIR_TUMOR)

            runnable: annotation_dir && !has_existing
            skip: true
                return meta
        }

    //
    // MODULE: LINX visualiser
    //
    // Create process input channel
    // channel: [ meta_linx, annotation_dir ]
    ch_linx_visualiser_inputs = ch_inputs_sorted.runnable
        .map { meta, annotation_dir ->

            def meta_linx = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
            ]

            return [meta_linx, annotation_dir]
        }

    // Run process
    LINX_VISUALISER(
        ch_linx_visualiser_inputs,
        genome_version,
        ensembl_data_resources,
    )

    ch_versions = ch_versions.mix(LINX_VISUALISER.out.versions)

    //
    // MODULE: gpgr LINX report
    //
    // Create process input channel
    // channel: [ meta_gpgr, annotation_dir, visualiser_dir ]
    ch_gpgr_linx_inputs = WorkflowOncoanalyser.groupByMeta(
        ch_inputs_sorted.runnable,
        WorkflowOncoanalyser.restoreMeta(LINX_VISUALISER.out.plots, ch_inputs),
    )
        .map { meta, annotation_dir, visualiser_dir ->

            def meta_gpgr_linx = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
            ]

            return [meta_gpgr_linx, annotation_dir, visualiser_dir]
        }

    // Run process
    LINXREPORT(
        ch_gpgr_linx_inputs,
    )

    ch_versions = ch_versions.mix(LINXREPORT.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, visualiser_dir ]
    ch_visualiser_dir_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(LINX_VISUALISER.out.plots, ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    visualiser_dir = ch_visualiser_dir_out // channel: [ meta, visualiser_dir ]

    versions       = ch_versions           // channel: [ versions.yml ]
}
