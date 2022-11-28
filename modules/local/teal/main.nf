process TEAL {
  //conda (params.enable_conda ? "bioconda::hmftools-teal=1.0.1" : null)
  container 'docker.io/scwatts/teal:1.0.1--2'

  input:
  tuple val(meta), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai), path(tumor_wgs_metrics), path(normal_wgs_metrics), path(cobalt_dir), path(purple_dir)

  output:
  tuple val(meta), path('teal/'), emit: teal_dir
  path 'versions.yml'           , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def tumor_args = normal_bam ? """
    -tumor ${meta.get(['sample_name', 'tumor'])}
    -tumor_bam ${tumor_bam}
    -tumor_wgs_metrics ${tumor_wgs_metrics}
  """ : ''
  def reference_args = normal_bam ? """
    -reference ${meta.get(['sample_name', 'normal'])}
    -reference_bam ${normal_bam}
    -reference_wgs_metrics ${normal_wgs_metrics}
  """ : ''

  """
  java \\
    -Xmx${task.memory.giga}g \\
    -jar "${task.ext.jarPath}" \\
      ${args} \\
      ${tumor_args.replaceAll('\n', '')} \\
      ${reference_args.replaceAll('\n', '')} \\
      -cobalt "${cobalt_dir}" \\
      -purple "${purple_dir}" \\
      -threads "${task.cpus}" \\
      -output_dir teal/

  # NOTE(SW): hard coded since there is no reliable way to obtain version information.
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      teal: 1.0.1
  END_VERSIONS
  """

  stub:
  """
  mkdir -p teal/
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
