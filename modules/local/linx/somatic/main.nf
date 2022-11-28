process LINX_SOMATIC {
  //conda (params.enable_conda ? "bioconda::hmftools-linx=1.21" : null)
  container 'docker.io/scwatts/linx:1.22--0'

  input:
  tuple val(meta), path(purple_dir)
  val genome_ver
  path fragile_sites
  path lines
  path ensembl_data_dir
  path known_fusion_data
  path driver_gene_panel

  output:
  tuple val(meta), path('linx_somatic/'), emit: annotation_dir
  path 'versions.yml'                   , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''

  """
  java \\
    -Xmx${task.memory.giga}g \\
    -jar "${task.ext.jarPath}" \\
      ${args} \\
      -sample "${meta.get(['sample_name', 'tumor'])}" \\
      -ref_genome_version "${genome_ver}" \\
      -sv_vcf "${purple_dir}/${meta.get(['sample_name', 'tumor'])}.purple.sv.vcf.gz" \\
      -purple_dir "${purple_dir}" \\
      -fragile_site_file "${fragile_sites}" \\
      -line_element_file "${lines}" \\
      -ensembl_data_dir "${ensembl_data_dir}" \\
      -check_fusions \\
      -known_fusion_file "${known_fusion_data}" \\
      -check_drivers \\
      -driver_gene_panel "${driver_gene_panel}" \\
      -output_dir linx_somatic/

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      linx: \$(java -jar "${task.ext.jarPath}" | sed 's/^.*LINX version: //')
  END_VERSIONS
  """

  stub:
  """
  mkdir linx_somatic/
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
