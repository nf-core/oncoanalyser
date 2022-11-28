//
// GRIPSS performs SV filtering.
//

include { GRIPSS_GERMLINE } from '../../modules/local/gripss/germline/main'
include { GRIPSS_SOMATIC  } from '../../modules/local/gripss/somatic/main'

workflow GRIPSS {
  take:
    ch_inputs               // channel: [val(meta), gridss_vcf]
    ref_data_genome_fa      //    file: /path/to/genome_fa
    ref_data_genome_fai     //    file: /path/to/genome_fai
    ref_data_genome_version //     val: genome version
    breakend_pon            //    file: /path/to/breakend_pon
    breakpoint_pon          //    file: /path/to/breakpoint_pon
    known_fusions           //    file: /path/to/known_fusions
    repeat_masker_file      //    file: /path/to/repeat_masker_file

  main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Germline
    GRIPSS_GERMLINE(
      ch_inputs,
      ref_data_genome_fa,
      ref_data_genome_fai,
      ref_data_genome_version,
      breakend_pon,
      breakpoint_pon,
      known_fusions,
      repeat_masker_file,
    )
    ch_versions = ch_versions.mix(GRIPSS_GERMLINE.out.versions)

    // Somatic
    GRIPSS_SOMATIC(
      ch_inputs,
      ref_data_genome_fa,
      ref_data_genome_fai,
      ref_data_genome_version,
      breakend_pon,
      breakpoint_pon,
      known_fusions,
      repeat_masker_file,
    )
    ch_versions = ch_versions.mix(GRIPSS_SOMATIC.out.versions)

  // Pack output
  // channel: [val(meta), hard_vcf, hard_tbi, soft_vcf, soft_tbi]
  ch_germline_out = WorkflowHmftools.group_by_meta(
    GRIPSS_GERMLINE.out.vcf_hard,
    GRIPSS_GERMLINE.out.vcf_soft,
  )
  // channel: [val(meta), hard_vcf, hard_tbi, soft_vcf, soft_tbi]
  ch_somatic_out = WorkflowHmftools.group_by_meta(
    GRIPSS_SOMATIC.out.vcf_hard,
    GRIPSS_SOMATIC.out.vcf_soft,
  )

  emit:
    germline = ch_germline_out // channel: [val(meta), hard_vcf, hard_tbi, soft_vcf, soft_tbi]
    somatic  = ch_somatic_out  // channel: [val(meta), hard_vcf, hard_tbi, soft_vcf, soft_tbi]

    versions = ch_versions     // channel: [versions.yml]
}
