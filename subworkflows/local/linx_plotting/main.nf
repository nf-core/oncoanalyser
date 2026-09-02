//
// LINX plotting visualises clusters structural variants
//

include { LINXREPORT      } from '../../../modules/local/linxreport/main'
include { LINX_VISUALISER } from '../../../modules/local/linx/visualiser/main'

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
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_linx_somatic_annotations,
        ch_amber_dir,
        ch_cobalt_dir,
        ch_purple_dir,
    )
        .map{ meta, linx_annotations_dir, amber_dir, cobalt_dir, purple_dir ->

            return [
                meta,
                Utils.selectCurrentOrExisting(linx_annotations_dir, meta, Constants.INPUT.LINX_ANNO_DIR_TUMOR),
                Utils.selectCurrentOrExisting(amber_dir, meta, Constants.INPUT.AMBER_DIR),
                Utils.selectCurrentOrExisting(cobalt_dir, meta, Constants.INPUT.COBALT_DIR),
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
            ]
        }
        .branch { meta, linx_annotations_dir, amber_dir, cobalt_dir, purple_dir ->

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.LINX_PLOT_DIR_TUMOR)

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
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
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
    ch_gpgr_linx_inputs = WorkflowOncoanalyser.groupByMeta(
        ch_inputs_sorted.runnable,
        WorkflowOncoanalyser.restoreMeta(channel.topic('linx_visualiser_plots'), ch_inputs),
    )
        .map { meta, linx_annotation_dir, _amber_dir, _cobalt_dir, _purple_dir, linx_visualiser_dir ->

            def meta_gpgr_linx = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
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
            WorkflowOncoanalyser.restoreMeta(channel.topic('linx_visualiser_plots'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    visualiser_dir = ch_outputs // channel: [ meta, linx_visualiser_dir ]
}
