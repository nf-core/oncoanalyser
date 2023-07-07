//
// LINX annotates and interprets structural variants
//

import Constants
import Utils

include { LINX_GERMLINE as GERMLINE } from '../../modules/local/linx/germline/main'
include { LINX_SOMATIC as SOMATIC   } from '../../modules/local/linx/somatic/main'

workflow LINX_ANNOTATION {
    take:
        // Sample data
        ch_inputs              // channel: [mandatory] [ meta ]
        ch_purple              // channel: [optional]  [ meta, purple_dir ]

        // Reference data
        genome_version         // channel: [mandatory] genome version
        ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/
        known_fusion_data      // channel: [mandatory] /path/to/known_fusion_data
        driver_gene_panel      // channel: [mandatory] /path/to/driver_gene_panel
        gene_id_file           // channel: [mandatory] /path/to/linx_gene_id_file

        // Params
        run_config             // channel: [mandatory] run configuration

    main:
        // Channel for versions.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [ meta, purple_dir ]
        ch_linx_inputs_source = run_config.stages.purple ? ch_purple : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PURPLE_DIR)

        //
        // MODULE: LINX germline annotation
        //
        // channel: [ meta, sv_vcf ]
        ch_linx_germline_out = Channel.empty()
        if (run_config.type == Constants.RunType.TUMOR_NORMAL) {

            // Create germline inputs and create process-specific meta
            // channel: [ meta_linx, sv_vcf ]
            ch_linx_inputs_germline = ch_linx_inputs_source
                .map { meta, purple_dir ->

                    def tumor_id = Utils.getTumorSampleName(meta, run_config.mode)
                    def sv_vcf = file(purple_dir).resolve("${tumor_id}.purple.sv.germline.vcf.gz")

                    if (!sv_vcf.exists()) {
                        return Constants.PLACEHOLDER_META
                    }

                    def meta_linx = [
                        key: meta.id,
                        id: Utils.getNormalWgsSampleName(meta),
                    ]

                    return [meta_linx, sv_vcf]
                }
                .filter { it != Constants.PLACEHOLDER_META }


            GERMLINE(
                ch_linx_inputs_germline,
                genome_version,
                ensembl_data_resources,
                driver_gene_panel,
            )

            // Set outputs, restoring original meta
            ch_versions = ch_versions.mix(GERMLINE.out.versions)
            ch_linx_germline_out = WorkflowOncoanalyser.restoreMeta(GERMLINE.out.annotation_dir, ch_inputs)
        }

        //
        // MODULE: LINX somatic annotation
        //
        // Create somatic inputs and create process-specific meta
        // channel: [ meta_linx, purple_dir ]
        ch_linx_inputs_somatic = ch_linx_inputs_source
            .map { meta, purple_dir ->
                def meta_linx = [
                    key: meta.id,
                    id: Utils.getTumorSampleName(meta, run_config.mode),
                ]
                return [meta_linx, purple_dir]
            }

        SOMATIC(
            ch_linx_inputs_somatic,
            genome_version,
            ensembl_data_resources,
            known_fusion_data,
            driver_gene_panel,
            gene_id_file,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(SOMATIC.out.versions)
        ch_linx_somatic_out = WorkflowOncoanalyser.restoreMeta(SOMATIC.out.annotation_dir, ch_inputs)

    emit:
        somatic       = ch_linx_somatic_out  // channel: [ meta, linx_annotation_dir ]
        germline      = ch_linx_germline_out // channel: [ meta, linx_annotation_dir ]

        versions      = ch_versions          // channel: [ versions.yml ]
}
