//
// LINX plotting visualises clusters structural variants
//

include { LINXREPORT      } from '../../../modules/local/linxreport/main'
include { LINX_VISUALISER } from '../../../modules/local/linx/visualiser/main'

workflow LINX_PLOTTING {
    take:
    // Sample data
    ch_inputs              // channel: [mandatory] [ meta ]
    ch_annotations         // channel: [mandatory] [ meta, annotation_dir ]

    // Sample data for copy number circos tracks
    ch_amber               // channel: [mandatory] [ meta, amber_dir ]
    ch_cobalt              // channel: [mandatory] [ meta, cobalt_dir ]
    ch_purple              // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_version         // channel: [mandatory] genome version
    ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/

    main:
    // Channel for versions.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources and sort
    // channel: runnable: [ meta, annotation_dir, amber_dir, cobalt_dir, purple_dir ]
    // channel: skip: [ meta ]

    ch_inputs_sorted = channels.WorkflowChannels.groupByMeta(
        ch_annotations,
        ch_amber,
        ch_cobalt,
        ch_purple,
    )
        .map{ meta, annotation_dir, amber_dir, cobalt_dir, purple_dir ->

            return [
                meta,
                Inputs.preferUserProvidedInput(annotation_dir, meta, sample.FileKey.LINX_ANNO_DIR_TUMOR),
                Inputs.preferUserProvidedInput(amber_dir, meta, sample.FileKey.AMBER_DIR),
                Inputs.preferUserProvidedInput(cobalt_dir, meta, sample.FileKey.COBALT_DIR),
                Inputs.preferUserProvidedInput(purple_dir, meta, sample.FileKey.PURPLE_DIR),
            ]
        }
        .branch { meta, annotation_dir, amber_dir, cobalt_dir, purple_dir ->

            def has_existing = Inputs.hasExistingInput(meta, sample.FileKey.LINX_PLOT_DIR_TUMOR)

            runnable: annotation_dir && !has_existing
            skip: true
                return meta
        }

    //
    // MODULE: LINX visualiser
    //
    // Create process input channel
    // channel: [ meta_linx, annotation_dir, amber_dir, cobalt_dir, purple_dir ]
    ch_linx_visualiser_inputs = ch_inputs_sorted.runnable
        .map { meta, annotation_dir, amber_dir, cobalt_dir, purple_dir ->

            def meta_linx = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Inputs.getTumorDnaSampleName(meta),
            ]

            return [meta_linx, annotation_dir, amber_dir, cobalt_dir, purple_dir]
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

    // NOTE(LN): UMCCR linx report code needs to be updated
    /*
    ch_gpgr_linx_inputs = channels.WorkflowChannels.groupByMeta(
        ch_inputs_sorted.runnable,
        channels.WorkflowChannels.restoreMeta(LINX_VISUALISER.out.plots, ch_inputs),
    )
        .map { meta, annotation_dir, amber_dir, cobalt_dir, purple_dir, visualiser_dir ->

            def meta_gpgr_linx = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Inputs.getTumorDnaSampleName(meta),
            ]

            return [meta_gpgr_linx, annotation_dir, visualiser_dir]
        }

    // Run process
    LINXREPORT(
        ch_gpgr_linx_inputs,
    )

    ch_versions = ch_versions.mix(LINXREPORT.out.versions)
    */

    // Set outputs, restoring original meta
    // channel: [ meta, visualiser_dir ]
    ch_visualiser_dir_out = Channel.empty()
        .mix(
            channels.WorkflowChannels.restoreMeta(LINX_VISUALISER.out.plots, ch_inputs),
            channels.PlaceholderChannels.toolDir(ch_inputs_sorted.skip),
        )

    emit:
    visualiser_dir = ch_visualiser_dir_out // channel: [ meta, visualiser_dir ]

    versions       = ch_versions           // channel: [ versions.yml ]
}
