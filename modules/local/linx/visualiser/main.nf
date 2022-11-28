process VISUALISER {
  //conda (params.enable_conda ? "bioconda::hmftools-linx=1.21" : null)
  container 'docker.io/scwatts/linx:1.22--0'

  // Lenient caching is required here as the hash of the LINX_SOMATIC
  // 'annotation_dir' output changes between runs even when the LINX_SOMATIC
  // process itself is cached.
  cache = 'lenient'

  input:
  tuple val(meta), path(linx)
  val genome_ver
  path ensembl_data_dir

  output:
  tuple val(meta), path('linx_visualiser/plot/'), emit: visualiser_dir
  path 'versions.yml'                           , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''

  """
  java \\
    -Xmx${task.memory.giga}g \\
    -cp "${task.ext.jarPath}" \\
    com.hartwig.hmftools.linx.visualiser.SvVisualiser \\
      ${args} \\
      -sample "${meta.get(['sample_name', 'tumor'])}" \\
      -ref_genome_version "${genome_ver}" \\
      -ensembl_data_dir "${ensembl_data_dir}" \\
      -plot_out linx_visualiser/plot \\
      -data_out linx_visualiser/data \\
      -vis_file_dir "${linx}" \\
      -circos "${task.ext.path_circos}" \\
      -threads "${task.cpus}"

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      linx: \$(java -jar "${task.ext.jarPath}" | sed 's/^.*LINX version: //')
  END_VERSIONS
  """

  stub:
  """
  mkdir -p linx_visualiser/plot/
  echo -e '${task.process}:\n  stub: noversions\n' > versions.yml
  """
}
