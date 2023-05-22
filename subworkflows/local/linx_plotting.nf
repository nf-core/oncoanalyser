//
// XXX
//
import Utils

include { GPGR_LINX as GPGR             } from '../../modules/local/gpgr/linx/main'
include { LINX_VISUALISER as VISUALISER } from '../../modules/local/linx/visualiser/main'

workflow LINX_PLOTTING {
    take:
        // Sample data
        ch_inputs
        ch_annotations

        // Reference data
        ref_data_genome_version         //     val: genome version
        ref_data_ensembl_data_resources //    file: /path/to/ensembl_data_resources/

        // Params
        run_config

    main:
        // Channel for versions.yml files
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [val(meta_linx), anno_dir]
        ch_linx_visualiser_inputs = ch_annotations
            .map { meta, anno_dir ->
                def meta_linx = [
                    key: meta.id,
                    id: Utils.getTumorSampleName(meta, run_config.mode),
                ]
                return [meta_linx, anno_dir]
            }

        VISUALISER(
            ch_linx_visualiser_inputs,
            ref_data_genome_version,
            ref_data_ensembl_data_resources,
        )

        // Set outputs, restoring original meta
        ch_visualiser_out = WorkflowOncoanalyser.restoreMeta(VISUALISER.out.visualiser_dir, ch_inputs)
        ch_versions = ch_versions.mix(VISUALISER.out.versions)


        // Create inputs and create process-specific meta
        // channel: [meta(meta_gpgr_linx), anno_dir, vis_dir]
        ch_gpgr_linx_inputs = WorkflowOncoanalyser.groupByMeta(
            ch_annotations,
            ch_visualiser_out,
        )
            .map { meta, anno_dir, vis_dir ->
                def meta_gpgr_linx = [
                    key: meta.id,
                    id: Utils.getTumorSampleName(meta, run_config.mode),
                ]
                return [meta_gpgr_linx, anno_dir, vis_dir]
            }

        GPGR(
            ch_gpgr_linx_inputs,
        )

        ch_versions = ch_versions.mix(GPGR.out.versions)

    emit:
        visualiser_dir = ch_visualiser_out // channel: [val(meta), visualiser_dir]

        versions = ch_versions             // channel: [versions.yml]
}
