//
// PAVE annotates somatic and germline variant VCFs with gene and transcript coding and protein effects.
//

include { PAVE_GERMLINE as GERMLINE } from '../../modules/local/pave/germline/main'
include { PAVE_SOMATIC as SOMATIC   } from '../../modules/local/pave/somatic/main'

workflow PAVE_ANNOTATION {
    take:

        // Sample data
        ch_inputs                       // channel: [val(meta)]
        ch_sage_germline_vcf            // channel: [val(meta), sage_germline_vcf]
        ch_sage_somatic_vcf             // channel: [val(meta), sage_somatic_vcf]

        // Reference data
        ref_data_genome_fasta           //    file: /path/to/genome_fasta
        ref_data_genome_fai             //    file: /path/to/genome_fai
        ref_data_genome_version         //     val: genome version
        ref_data_sage_pon               //    file: /path/to/sage_pon
        ref_data_sage_blocklist_regions //    file: /path/to/sage_blocklist_regions
        ref_data_sage_blocklist_sites   //    file: /path/to/sage_blocklist_sites
        ref_data_clinvar_annotations    //    file: /path/to/clinvar_annotations
        ref_data_segment_mappability    //    file: /path/to/segment_mappability
        ref_data_driver_gene_panel      //    file: /path/to/driver_gene_panel
        ref_data_ensembl_data_resources //    file: /path/to/ensembl_data_resources/
        ref_data_gnomad_pon_dir         //    file: /path/to/gnomad_pon_dir/

        // Params
        run

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [meta, sage_vcf, sage_tbi]
        if (run.sage) {
            ch_pave_germline_inputs_source = ch_sage_germline_vcf
            ch_pave_somatic_inputs_source = ch_sage_somatic_vcf
        } else {
            ch_pave_germline_inputs_source = WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.SAGE_VCF_TUMOR)
                .filter { it[0] != Constants.META_PLACEHOLDER }
            ch_pave_somatic_inputs_source = WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.SAGE_VCF_NORMAL)
                .filter { it[0] != Constants.META_PLACEHOLDER }
        }

        // Create inputs and create process-specific meta
        // channel: [val(meta_pave), sage_vcf]
        ch_pave_germline_inputs = ch_pave_germline_inputs_source
            .map { meta, sage_vcf, sage_tbi ->
                def pave_meta = [
                    key: meta.id,
                    // NOTE(SW): use of tumor sample name for PAVE germline is correct
                    id: Utils.getTumorWgsSampleName(meta),
                ]
                return [pave_meta, sage_vcf]
            }

        ch_pave_somatic_inputs = ch_pave_somatic_inputs_source
            .map { meta, sage_vcf, sage_tbi ->
                def pave_meta = [
                    key: meta.id,
                    id: Utils.getTumorWgsSampleName(meta),
                ]
                return [pave_meta, sage_vcf]
            }

        GERMLINE(
            ch_pave_germline_inputs,
            ref_data_genome_fasta,
            ref_data_genome_fai,
            ref_data_genome_version,
            ref_data_sage_blocklist_regions,
            ref_data_sage_blocklist_sites,
            ref_data_clinvar_annotations,
            ref_data_segment_mappability,
            ref_data_driver_gene_panel,
            ref_data_ensembl_data_resources,
        )

        SOMATIC(
            ch_pave_somatic_inputs,
            ref_data_genome_fasta,
            ref_data_genome_fai,
            ref_data_genome_version,
            ref_data_sage_pon,
            ref_data_segment_mappability,
            ref_data_driver_gene_panel,
            ref_data_ensembl_data_resources,
            ref_data_gnomad_pon_dir,
        )

        // Set outputs
        ch_germline_out = WorkflowOncoanalyser.restoreMeta(GERMLINE.out.vcf, ch_inputs)
        ch_somatic_out = WorkflowOncoanalyser.restoreMeta(SOMATIC.out.vcf, ch_inputs)
        ch_versions = ch_versions.mix(
            GERMLINE.out.versions,
            SOMATIC.out.versions,
        )

    emit:
        germline = ch_germline_out // channel: [val(meta), pave_vcf]
        somatic = ch_somatic_out   // channel: [val(meta), pave_vcf]

        versions = ch_versions           // channel: [versions.yml]
}
