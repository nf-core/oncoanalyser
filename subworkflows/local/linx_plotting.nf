//
// LINX plotting visualises clusters structural variants
//

import Constants
import Utils

include { GPGR_LINX as GPGR             } from '../../modules/local/gpgr/linx/main'
include { LINX_VISUALISER as VISUALISER } from '../../modules/local/linx/visualiser/main'

workflow LINX_PLOTTING {
    take:
        // Sample data
        ch_inputs                     // channel: [mandatory] [ meta ]
        ch_annotations                // channel: [mandatory] [ meta, annotation_dir ]

        // Reference data
        genome_version                // channel: [mandatory] genome version
        ensembl_data_resources        // channel: [mandatory] /path/to/ensembl_data_resources/

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
        VISUALISER(
            ch_linx_visualiser_inputs,
            genome_version,
            ensembl_data_resources,
        )

        ch_versions = ch_versions.mix(VISUALISER.out.versions)

        //
        // MODULE: gpgr LINX report
        //
        // Create process input channel
        // channel: [ meta_gpgr, annotation_dir, visualiser_dir_all ]
        ch_gpgr_linx_inputs = WorkflowOncoanalyser.groupByMeta(
            ch_inputs_sorted.runnable,
            WorkflowOncoanalyser.restoreMeta(VISUALISER.out.visualiser_dir_all, ch_inputs),
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
        GPGR(
            ch_gpgr_linx_inputs,
        )

        ch_versions = ch_versions.mix(GPGR.out.versions)

        // Set outputs, restoring original meta
        // channel: [ meta, visualiser_dir_all ]
        ch_visualiser_all_out = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(VISUALISER.out.visualiser_dir_all, ch_inputs),
                ch_inputs_sorted.skip.map { meta -> [meta, []] },
            )

        // channel: [ meta, visualiser_dir_reportable ]
        ch_visualiser_reportable_out = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(VISUALISER.out.visualiser_dir_reportable, ch_inputs),
                ch_inputs_sorted.skip.map { meta -> [meta, []] },
            )

    emit:
        visualiser_dir_all        = ch_visualiser_all_out        // channel: [ meta, visualiser_dir_all ]
        visualiser_dir_reportable = ch_visualiser_reportable_out // channel: [ meta, visualiser_dir_reportable ]

        versions = ch_versions                                   // channel: [ versions.yml ]
}
