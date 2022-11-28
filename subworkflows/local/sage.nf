//
// SAGE is a precise and highly sensitive somatic SNV, MNV and small INDEL caller.
//

include { SAGE_GERMLINE } from '../../modules/local/sage/germline/main'
include { SAGE_SOMATIC  } from '../../modules/local/sage/somatic/main'

workflow SAGE {
  take:
    ch_inputs                             // channel: [meta, tumor_bam, normal_bam, tumor_bai, normal_bai]
    ref_data_genome_fa                    //    file: /path/to/genome_fa
    ref_data_genome_fai                   //    file: /path/to/genome_fai
    ref_data_genome_dict                  //    file: /path/to/genome_dict
    ref_data_genome_version               //     val: genome version
    ref_data_sage_known_hotspots_germline //    file: /path/to/sage_known_hotspots_germline
    ref_data_sage_known_hotspots_somatic  //    file: /path/to/sage_known_hotspots_somatic
    ref_data_sage_coding_panel            //    file: /path/to/sage_coding_panel
    ref_data_sage_high_confidence         //    file: /path/to/sage_high_confidence
    ref_data_sage_pon_file                //    file: /path/to/sage_pon
    ref_data_mappability_bed              //    file: /path/to/mappability_bed
    ref_data_driver_gene_panel            //    file: /path/to/driver_gene_panel
    ref_data_ensembl_data_dir             //    file: /path/to/ensembl_data_dir/

  main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Germline
    SAGE_GERMLINE(
      ch_inputs,
      ref_data_genome_fa,
      ref_data_genome_fai,
      ref_data_genome_dict,
      ref_data_genome_version,
      ref_data_sage_known_hotspots_germline,
      ref_data_sage_coding_panel,
      ref_data_sage_high_confidence,
      ref_data_ensembl_data_dir,
    )
    ch_versions = ch_versions.mix(SAGE_GERMLINE.out.versions)

    // Somatic
    SAGE_SOMATIC(
      ch_inputs,
      ref_data_genome_fa,
      ref_data_genome_fai,
      ref_data_genome_dict,
      ref_data_genome_version,
      ref_data_sage_known_hotspots_somatic,
      ref_data_sage_coding_panel,
      ref_data_sage_high_confidence,
      ref_data_ensembl_data_dir,
    )
    ch_versions = ch_versions.mix(SAGE_SOMATIC.out.versions)

  emit:
    germline = SAGE_GERMLINE.out.vcf // channel: [val(meta), sage_vcf]
    somatic = SAGE_SOMATIC.out.vcf   // channel: [val(meta), sage_vcf]

    versions = ch_versions           // channel: [versions.yml]
}
