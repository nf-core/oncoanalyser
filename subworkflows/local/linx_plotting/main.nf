//
// LINX plotting visualises clusters structural variants
//

nextflow.enable.types = true

include { LINXREPORT       } from '../../../modules/local/linxreport/main'
include { LINX_VISUALISER  } from '../../../modules/local/linx/visualiser/main'

include { getInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSample        } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSampleName    } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { groupByMeta              } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { joinMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { restoreMeta              } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { FileType                 } from '../utils_nfcore_oncoanalyser_pipeline/types_enums'
include { selectCurrentOrExisting  } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow LINX_PLOTTING {
    take:
    // Sample data
    ch_inputs                  : Channel<Map>              // channel: [mandatory] [ meta ]
    ch_linx_somatic_annotations: Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, linx_annotation_dir ]
    ch_amber_dir               : Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, amber_dir ]
    ch_cobalt_dir              : Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, cobalt_dir ]
    ch_purple_dir              : Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_version             : Channel<String>           // channel: [mandatory] genome version
    ensembl_data_resources     : Channel<Path>             // channel: [mandatory] /path/to/ensembl_data_resources/

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

            def meta_linx = record(
                key: meta.case_id,
                id: meta.case_id,
                sample_id: getTumorDnaSampleName(meta),
            )

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
    // NOTE(SW): subscribe once and reuse, channel.topic() is single-consumer under types
    ch_linx_plots_topic = channel.topic('linx_visualiser_plots')

    // channel: [ meta_gpgr, linx_annotation_dir, linx_visualiser_dir ]
    ch_gpgr_linx_inputs = groupByMeta([
        ch_inputs_sorted.runnable,
        restoreMeta(ch_linx_plots_topic, ch_inputs),
    ])
        .map { meta, linx_annotation_dir, _amber_dir, _cobalt_dir, _purple_dir, linx_visualiser_dir ->

            def meta_gpgr_linx = record(
                key: meta.case_id,
                id: meta.case_id,
                sample_id: getTumorDnaSampleName(meta),
            )

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
    ch_outputs_visualiser_dir = channel.empty()
        .mix(
            restoreMeta(ch_linx_plots_topic, ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    // channel: [ meta, linxreport_html ]
    ch_outputs_linxreport_html = channel.empty()
        .mix(
            restoreMeta(channel.topic('linxreport_html'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    emit:
    visualiser_dir  = ch_outputs_visualiser_dir  // channel: [ meta, linx_visualiser_dir ]
    linxreport_html = ch_outputs_linxreport_html // channel: [ meta, linxreport_html ]
}
