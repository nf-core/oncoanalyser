//
// Linx is an annotation, interpretation and visualisation tool for structural variants.
//
import Constants
import Utils

include { LINX_GERMLINE as GERMLINE } from '../../modules/local/linx/germline/main'
include { LINX_SOMATIC as SOMATIC   } from '../../modules/local/linx/somatic/main'

workflow LINX_ANNOTATION {
    take:
        // Sample data
        ch_inputs                       // channel: [val(meta)]
        ch_purple                       // channel: [val(meta), purple_dir]

        // Reference data
        ref_data_genome_version         //     val: genome version
        ref_data_ensembl_data_resources //    file: /path/to/ensembl_data_resources/
        ref_data_known_fusion_data      //    file: /path/to/known_fusion_data
        ref_data_driver_gene_panel      //    file: /path/to/driver_gene_panel
        gene_id_file                    //    file: /path/to/linx_gene_id_file

        // Params
        run

    main:
        // Channel for versions.yml files
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [val(meta), purple_dir]
        ch_linx_inputs_source = run.purple ? ch_purple : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PURPLE_DIR)

        // Create inputs and create process-specific meta
        // channel: [val(meta_linx), sv_vcf]
        ch_linx_inputs_germline = ch_linx_inputs_source
            .map { meta, purple_dir ->

                def tumor_id = Utils.getTumorWgsSampleName(meta)
                def sv_vcf = file(purple_dir).resolve("${tumor_id}.purple.sv.vcf.gz")

                if (!sv_vcf.exists()) {
                    return Constants.META_PLACEHOLDER
                }

                def meta_linx = [
                    key: meta.id,
                    id: Utils.getNormalWgsSampleName(meta),
                ]

                return [meta_linx, sv_vcf]
            }
            .filter { it != Constants.META_PLACEHOLDER }

        // channel: [val(meta_linx), purple_dir]
        ch_linx_inputs_somatic = ch_linx_inputs_source
            .map { meta, purple_dir ->
                def meta_linx = [
                    key: meta.id,
                    id: Utils.getTumorWgsSampleName(meta),
                ]
                return [meta_linx, purple_dir]
            }

        GERMLINE(
            ch_linx_inputs_germline,
            ref_data_genome_version,
            ref_data_ensembl_data_resources,
            ref_data_driver_gene_panel,
        )

        SOMATIC(
            ch_linx_inputs_somatic,
            ref_data_genome_version,
            ref_data_ensembl_data_resources,
            ref_data_known_fusion_data,
            ref_data_driver_gene_panel,
            gene_id_file,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(
            SOMATIC.out.versions,
            GERMLINE.out.versions,
        )

        // channel: [val(meta), linx_annotation_dir]
        ch_linx_somatic_out = WorkflowOncoanalyser.restoreMeta(SOMATIC.out.annotation_dir, ch_inputs)
        ch_linx_germline_out = WorkflowOncoanalyser.restoreMeta(GERMLINE.out.annotation_dir, ch_inputs)

    emit:
        somatic       = ch_linx_somatic_out       // channel: [val(meta), linx_annotation_dir]
        germline      = ch_linx_germline_out      // channel: [val(meta), linx_annotation_dir]

        versions      = ch_versions               // channel: [versions.yml]
}
