//
// PAVE annotates somatic and germline variant VCFs with gene and transcript coding and protein effects
//

include { PAVE_GERMLINE as GERMLINE } from '../../modules/local/pave/germline/main'
include { PAVE_SOMATIC as SOMATIC   } from '../../modules/local/pave/somatic/main'

workflow PAVE_ANNOTATION {
    take:
        // Sample data
        ch_inputs              // channel: [mandatory] [ meta ]
        ch_sage_germline_vcf   // channel: [optional]  [ meta, sage_germline_vcf ]
        ch_sage_somatic_vcf    // channel: [mandatory] [ meta, sage_somatic_vcf ]

        // Reference data
        genome_fasta           // channel: [mandatory] /path/to/genome_fasta
        genome_version         // channel: [mandatory] genome version
        genome_fai             // channel: [mandatory] /path/to/genome_fai
        sage_pon               // channel: [mandatory] /path/to/sage_pon
        sage_pon_artefacts     // channel: [mandatory] /path/to/sage_pon_artefacts
        sage_blocklist_regions // channel: [mandatory] /path/to/sage_blocklist_regions
        sage_blocklist_sites   // channel: [mandatory] /path/to/sage_blocklist_sites
        clinvar_annotations    // channel: [mandatory] /path/to/clinvar_annotations
        segment_mappability    // channel: [mandatory] /path/to/segment_mappability
        driver_gene_panel      // channel: [mandatory] /path/to/driver_gene_panel
        ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/
        gnomad_resource        // channel: [mandatory] /path/to/gnomad_resource

        // Params
        run_config             // channel: [mandatory] run configuration

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [meta, sage_vcf, sage_tbi]
        if (run_config.stages.sage) {
            ch_pave_germline_inputs_source = ch_sage_germline_vcf
            ch_pave_somatic_inputs_source = ch_sage_somatic_vcf
        } else {
            ch_pave_germline_inputs_source = WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.SAGE_VCF_NORMAL)
                .filter { it[0] != Constants.PLACEHOLDER_META }
                .map { meta, vcf -> [meta, vcf, "${vcf}.tbi"] }
            ch_pave_somatic_inputs_source = WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.SAGE_VCF_TUMOR)
                .filter { it[0] != Constants.PLACEHOLDER_META }
                .map { meta, vcf -> [meta, vcf, "${vcf}.tbi"] }
        }

        //
        // MODULE: PAVE germline
        //
        ch_germline_out = Channel.empty()
        if (run_config.type == Constants.RunType.TUMOR_NORMAL) {

            // Create inputs and create process-specific meta
            // channel: [val(meta_pave), sage_vcf]
            ch_pave_germline_inputs = ch_pave_germline_inputs_source
                .map { meta, sage_vcf, sage_tbi ->
                    def pave_meta = [
                        key: meta.id,
                        // NOTE(SW): use of tumor sample name for PAVE germline is correct
                        id: Utils.getTumorSampleName(meta, run_config.mode),
                    ]
                    return [pave_meta, sage_vcf]
                }

            GERMLINE(
                ch_pave_germline_inputs,
                genome_fasta,
                genome_version,
                genome_fai,
                sage_blocklist_regions,
                sage_blocklist_sites,
                clinvar_annotations,
                segment_mappability,
                driver_gene_panel,
                ensembl_data_resources,
            )

            ch_versions = ch_versions.mix(GERMLINE.out.versions)
            ch_germline_out = WorkflowOncoanalyser.restoreMeta(GERMLINE.out.vcf, ch_inputs)

        }

        //
        // MODULE: PAVE somatic
        //
        // Create inputs and create process-specific meta
        // channel: [val(meta_pave), sage_vcf]
        ch_pave_somatic_inputs = ch_pave_somatic_inputs_source
            .map { meta, sage_vcf, sage_tbi ->
                def pave_meta = [
                    key: meta.id,
                    id: Utils.getTumorSampleName(meta, run_config.mode),
                ]
                return [pave_meta, sage_vcf]
            }

        SOMATIC(
            ch_pave_somatic_inputs,
            genome_fasta,
            genome_version,
            genome_fai,
            sage_pon,
            sage_pon_artefacts,
            segment_mappability,
            driver_gene_panel,
            ensembl_data_resources,
            gnomad_resource,
        )

        ch_versions = ch_versions.mix(SOMATIC.out.versions)
        ch_somatic_out = WorkflowOncoanalyser.restoreMeta(SOMATIC.out.vcf, ch_inputs)

    emit:
        germline = ch_germline_out // channel: [val(meta), pave_vcf]
        somatic = ch_somatic_out   // channel: [val(meta), pave_vcf]

        versions = ch_versions     // channel: [versions.yml]
}
