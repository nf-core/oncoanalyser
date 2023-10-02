//
// LINX plotting visualises clusters structural variants
//

import Utils

include { GPGR_LINX as GPGR             } from '../../modules/local/gpgr/linx/main'
include { LINX_VISUALISER as VISUALISER } from '../../modules/local/linx/visualiser/main'

workflow LINX_PLOTTING {
    take:
        // Sample data
        ch_inputs              // channel: [mandatory] [ meta ]
        ch_annotations         // channel: [mandatory] [ meta, linx_annotation_dir ]

        // Reference data
        genome_version         // channel: [mandatory] genome version
        ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/

        // Params
        run_config             // channel: [mandatory] run configuration

    main:
        // Channel for versions.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [ meta_linx, linx_annotation_dir ]
        ch_linx_visualiser_inputs = ch_annotations
            .map { meta, anno_dir ->
                def meta_linx = [
                    key: meta.id,
                    id: Utils.getTumorDnaSampleName(meta),
                ]
                return [meta_linx, anno_dir]
            }

        VISUALISER(
            ch_linx_visualiser_inputs,
            genome_version,
            ensembl_data_resources,
        )

        // Set outputs, restoring original meta
        ch_visualiser_out = WorkflowOncoanalyser.restoreMeta(VISUALISER.out.visualiser_dir, ch_inputs)
        ch_versions = ch_versions.mix(VISUALISER.out.versions)


        // Create inputs and create process-specific meta
        // channel: [ meta_gpgr_linx, linx_annotation_dir, visualiser_dir ]
        ch_gpgr_linx_inputs = WorkflowOncoanalyser.groupByMeta(
            ch_annotations,
            ch_visualiser_out,
        )
            .map { meta, anno_dir, vis_dir ->
                def meta_gpgr_linx = [
                    key: meta.id,
                    id: Utils.getTumorDnaSampleName(meta),
                ]
                return [meta_gpgr_linx, anno_dir, vis_dir]
            }

        GPGR(
            ch_gpgr_linx_inputs,
        )

        ch_versions = ch_versions.mix(GPGR.out.versions)

    emit:
        visualiser_dir = ch_visualiser_out // channel: [ meta, visualiser_dir ]

        versions = ch_versions             // channel: [ versions.yml ]
}
