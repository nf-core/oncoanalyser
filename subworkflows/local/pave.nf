//
// PAVE annotates somatic and germline variant VCFs with gene and transcript coding and protein effects.
//

include { PAVE_GERMLINE } from '../../modules/local/pave/germline/main'
include { PAVE_SOMATIC  } from '../../modules/local/pave/somatic/main'

workflow PAVE {
  take:
    ch_inputs_germline          // channel: [val(meta), sage_germline_vcf]
    ch_inputs_somatic           // channel: [val(meta), sage_somatic_vcf]
    ref_data_genome_fa          //    file: /path/to/genome_fa
    ref_data_genome_fai         //    file: /path/to/genome_fai
    ref_data_genome_version     //     val: genome version
    ref_data_sage_pon_file      //    file: /path/to/sage_pon_file
    ref_data_sage_blacklist_bed //    file: /path/to/sage_blacklist_bed
    ref_data_sage_blacklist_vcf //    file: /path/to/sage_black_list_vcf
    ref_data_clinvar_vcf        //    file: /path/to/clinvar_vcf
    ref_data_mappability_bed    //    file: /path/to/mappability_bed
    ref_data_driver_gene_panel  //    file: /path/to/driver_gene_panel
    ref_data_ensembl_data_dir   //    file: /path/to/ensembl_data_dir/

  main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Germline
    PAVE_GERMLINE(
      ch_inputs_germline,
      ref_data_genome_fa,
      ref_data_genome_fai,
      ref_data_genome_version,
      ref_data_sage_blacklist_bed,
      ref_data_sage_blacklist_vcf,
      ref_data_clinvar_vcf,
      ref_data_mappability_bed,
      ref_data_driver_gene_panel,
      ref_data_ensembl_data_dir,
    )
    ch_versions = ch_versions.mix(PAVE_GERMLINE.out.versions)

    // Somatic
    PAVE_SOMATIC(
      ch_inputs_somatic,
      ref_data_genome_fa,
      ref_data_genome_fai,
      ref_data_genome_version,
      ref_data_sage_pon_file,
      ref_data_mappability_bed,
      ref_data_driver_gene_panel,
      ref_data_ensembl_data_dir,
    )
    ch_versions = ch_versions.mix(PAVE_SOMATIC.out.versions)

  emit:
    germline = PAVE_GERMLINE.out.vcf // channel: [val(meta), pave_vcf]
    somatic = PAVE_SOMATIC.out.vcf   // channel: [val(meta), pave_vcf]

    versions = ch_versions           // channel: [versions.yml]
}
