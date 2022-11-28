process AMBER {
  //conda (params.enable_conda ? "bioconda::hmftools-amber=3.9" : null)
  container 'docker.io/scwatts/amber:3.9--3'

  input:
  tuple val(meta), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai)
  val ref_genome_ver
  path loci

  output:
  tuple val(meta), path('amber/'), emit: amber_dir
  path 'versions.yml'            , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''

  """
  java \\
    -Xmx${task.memory.giga}g \\
    -jar "${task.ext.jarPath}" \\
      ${args} \\
      -tumor "${meta.get(['sample_name', 'tumor'])}" \\
      -tumor_bam "${tumor_bam}" \\
      -reference "${meta.get(['sample_name', 'normal'])}" \\
      -reference_bam "${normal_bam}" \\
      -ref_genome_version ${ref_genome_ver} \\
      -output_dir amber/ \\
      -threads "${task.cpus}" \\
      -loci "${loci}"

  # NOTE(SW): hard coded since there is no reliable way to obtain version information.
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      amber: 3.9
  END_VERSIONS
  """

  stub:
  """
  mkdir -p amber/
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
