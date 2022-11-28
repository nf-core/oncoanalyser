//
// Linx is an annotation, interpretation and visualisation tool for structural variants.
//

include { LINX_GERMLINE } from '../../modules/local/linx/germline/main'
include { LINX_SOMATIC  } from '../../modules/local/linx/somatic/main'
include { VISUALISER    } from '../../modules/local/linx/visualiser/main'

workflow LINX {
  take:
    ch_linx_germline_inputs     // channel: [val(meta), gripss_hard_vcf]
    ch_linx_somatic_inputs      // channel: [val(meta), purple_dir]
    ref_data_genome_version     //     val: genome version
    ref_data_linx_fragile_sites //    file: /path/to/linx_fragile_sites
    ref_data_linx_lines         //    file: /path/to/linx_lines
    ref_data_ensembl_data_dir   //    file: /path/to/ensembl_data_dir/
    ref_data_known_fusion_data  //    file: /path/to/known_fusion_data
    ref_data_driver_gene_panel  //    file: /path/to/driver_gene_panel

  main:
    // Channel for versions.yml files
    ch_versions = Channel.empty()

    // Germline
    LINX_GERMLINE(
      ch_linx_germline_inputs,
      ref_data_genome_version,
      ref_data_linx_fragile_sites,
      ref_data_linx_lines,
      ref_data_ensembl_data_dir,
      ref_data_driver_gene_panel,
    )
    ch_versions = ch_versions.mix(LINX_GERMLINE.out.versions)

    // Somatic
    LINX_SOMATIC(
      ch_linx_somatic_inputs,
      ref_data_genome_version,
      ref_data_linx_fragile_sites,
      ref_data_linx_lines,
      ref_data_ensembl_data_dir,
      ref_data_known_fusion_data,
      ref_data_driver_gene_panel,
    )
    ch_versions = ch_versions.mix(LINX_SOMATIC.out.versions)

    VISUALISER(
      LINX_SOMATIC.out.annotation_dir,
      ref_data_genome_version,
      ref_data_ensembl_data_dir,
    )
    ch_versions = ch_versions.mix(VISUALISER.out.versions)

    // channel: [val(meta), linx_annotation_dir, linx_visualiser_dir]
    ch_linx_somatic_out = WorkflowHmftools.group_by_meta(
      LINX_SOMATIC.out.annotation_dir,
      VISUALISER.out.visualiser_dir,
    )

  emit:
    somatic  = ch_linx_somatic_out              // channel: [val(meta), linx_annotation_dir, linx_visualiser_dir]
    germline = LINX_GERMLINE.out.annotation_dir // channel: [val(meta), linx_annotation_dir]

    versions = ch_versions                      // channel: [versions.yml]
}
