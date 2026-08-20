//
// LINX plotting visualises clusters structural variants
//

include { LINXREPORT  } from '../../../modules/local/linxreport/main'
include { LINX_VISUALISER  } from '../../../modules/local/linx/visualiser/main'

include { FileType                 } from '../utils_nfcore_oncoanalyser_pipeline/types'
include { groupByMeta              } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { joinMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { restoreMeta              } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { getInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSample        } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSampleName    } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { hasInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { selectCurrentOrExisting  } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow LINX_PLOTTING {
    take:
    // Sample data
    ch_inputs                   // channel: [mandatory] [ meta ]
    ch_linx_somatic_annotations // channel: [mandatory] [ meta, linx_annotation_dir ]
    ch_amber_dir                // channel: [mandatory] [ meta, amber_dir ]
    ch_cobalt_dir               // channel: [mandatory] [ meta, cobalt_dir ]
    ch_purple_dir               // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_version              // channel: [mandatory] genome version
    ensembl_data_resources      // channel: [mandatory] /path/to/ensembl_data_resources/

    main:
    //
    // STEP: Handle inputs
    //
    // Select input sources then sort
    // channel: runnable: [ meta, linx_annotation_dir, amber_dir, cobalt_dir, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = groupByMeta([
        ch_linx_somatic_annotations,
        ch_amber_dir,
        ch_cobalt_dir,
        ch_purple_dir,
    ])
        .map{ meta, linx_annotations_dir, amber_dir, cobalt_dir, purple_dir ->

            return [
                meta,
                selectCurrentOrExisting(linx_annotations_dir, getInput(getTumorDnaSample(meta), FileType.LINX_ANNO_DIR)),
                selectCurrentOrExisting(amber_dir, getInput(getTumorDnaSample(meta), FileType.AMBER_DIR)),
                selectCurrentOrExisting(cobalt_dir, getInput(getTumorDnaSample(meta), FileType.COBALT_DIR)),
                selectCurrentOrExisting(purple_dir, getInput(getTumorDnaSample(meta), FileType.PURPLE_DIR)),
            ]
        }
        .branch { meta, linx_annotations_dir, amber_dir, cobalt_dir, purple_dir ->

            def has_existing = hasInput(getTumorDnaSample(meta), FileType.LINX_PLOT_DIR)

            runnable: linx_annotations_dir && ! has_existing
            skip: true
                return meta
        }

    //
    // MODULE: LINX visualiser
    //
    // Create process input channel
    // channel: [ meta_linx, linx_annotation_dir, amber_dir, cobalt_dir, purple_dir ]
    ch_linx_visualiser_inputs = ch_inputs_sorted.runnable
        .map { meta, linx_annotations_dir, amber_dir, cobalt_dir, purple_dir ->

            def meta_linx = [
                key: meta.case_id,
                id: meta.case_id,
                sample_id: getTumorDnaSampleName(meta),
            ]

            return [meta_linx, linx_annotations_dir, amber_dir, cobalt_dir, purple_dir]
        }

    // Run process
    LINX_VISUALISER(
        ch_linx_visualiser_inputs,
        genome_version,
        ensembl_data_resources,
    )

    //
    // MODULE: gpgr LINX report
    //
    // Create process input channel
    // channel: [ meta_gpgr, linx_annotation_dir, linx_visualiser_dir ]
    ch_gpgr_linx_inputs = groupByMeta([
        ch_inputs_sorted.runnable,
        restoreMeta(channel.topic('linx_visualiser_plots'), ch_inputs),
    ])
        .map { meta, linx_annotation_dir, _amber_dir, _cobalt_dir, _purple_dir, linx_visualiser_dir ->

            def meta_gpgr_linx = [
                key: meta.case_id,
                id: meta.case_id,
                sample_id: getTumorDnaSampleName(meta),
            ]

            return [meta_gpgr_linx, linx_annotation_dir, linx_visualiser_dir]
        }

    // Run process
    LINXREPORT(
        ch_gpgr_linx_inputs,
    )

    //
    // STEP: Handle outputs
    //
    // Set outputs, restoring original meta
    // channel: [ meta, linx_visualiser_dir ]
    ch_outputs = channel.empty()
        .mix(
            restoreMeta(channel.topic('linx_visualiser_plots'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    visualiser_dir = ch_outputs // channel: [ meta, linx_visualiser_dir ]
}
