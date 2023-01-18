//
// PAVE annotates somatic and germline variant VCFs with gene and transcript coding and protein effects.
//

include { PAVE_GERMLINE } from '../../modules/local/pave/germline/main'
include { PAVE_SOMATIC  } from '../../modules/local/pave/somatic/main'

workflow PAVE {
    take:
        ch_inputs_germline              // channel: [val(meta), sage_germline_vcf]
        ch_inputs_somatic               // channel: [val(meta), sage_somatic_vcf]
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

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Germline
        PAVE_GERMLINE(
            ch_inputs_germline,
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
        ch_versions = ch_versions.mix(PAVE_GERMLINE.out.versions)

        // Somatic
        PAVE_SOMATIC(
            ch_inputs_somatic,
            ref_data_genome_fasta,
            ref_data_genome_fai,
            ref_data_genome_version,
            ref_data_sage_pon,
            ref_data_segment_mappability,
            ref_data_driver_gene_panel,
            ref_data_ensembl_data_resources,
            ref_data_gnomad_pon_dir,
        )
        ch_versions = ch_versions.mix(PAVE_SOMATIC.out.versions)

    emit:
        germline = PAVE_GERMLINE.out.vcf // channel: [val(meta), pave_vcf]
        somatic = PAVE_SOMATIC.out.vcf   // channel: [val(meta), pave_vcf]

        versions = ch_versions           // channel: [versions.yml]
}
