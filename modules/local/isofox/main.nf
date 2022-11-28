// NOTE(SW): care must be taken for using exp_counts; only supports read length of 151 bp


process ISOFOX {
  container 'docker.io/scwatts/isofox:1.5--0'

  input:
  tuple val(meta), path(bam)
  path genome_fa
  path genome_fai
  val genome_ver
  path ensembl_data_dir
  path exp_counts
  path exp_gc_ratios

  output:
  tuple val(meta), path('isofox/'), emit: isofox_dir
  path 'versions.yml'             , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''

  """
  mkdir -p isofox/

  java \\
    -Xmx${task.memory.giga}g \\
    -jar ${task.ext.jarPath} \\
      ${args} \\
      -sample ${meta.get(['sample_name', 'tumor'])} \\
      -bam_file ${bam} \\
      -ref_genome ${genome_fa} \\
      -ref_genome_version ${genome_ver} \\
      -ensembl_data_dir ${ensembl_data_dir} \\
      -exp_counts_file ${exp_counts} \\
      -exp_gc_ratios_file ${exp_gc_ratios} \\
      -output_dir ./isofox/ \\
      -threads ${task.cpus}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      isofox: \$(java -jar ${task.ext.jarPath} | sed -n '1s/^.*version: //p')
  END_VERSIONS
  """

  stub:
  """
  mkdir -p isofox/
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
