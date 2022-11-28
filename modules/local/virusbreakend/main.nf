// NOTE(SW): the --db argument for the virusbreakend command must have a trailing slash if it is a
// symlink

process VIRUSBREAKEND {
  container 'docker.io/scwatts/gridss:2.13.2--3'

  input:
  tuple val(meta), path(bam)
  path virusbreakenddb
  path gridss_config
  path genome_fa
  path genome_fai
  path genome_dict
  path genome_bwa_index_dir, stageAs: 'bwa_index'
  path genome_bwa_index_image
  path genome_gridss_index

  output:
  tuple val(meta), path("*.summary.tsv"), emit: tsv
  path "*.virusbreakend.vcf"            , emit: vcf
  path 'versions.yml'                   , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''

  """
  virusbreakend \\
    --jar ${task.ext.jarPath} \\
    --gridssargs "--jvmheap ${task.memory.giga}g" \\
    --threads ${task.cpus} \\
    --db ${virusbreakenddb.toString().replaceAll("/\$", "")}/ \\
    --output ${meta.id}.virusbreakend.vcf \\
    --reference ${genome_fa} \\
    ${bam}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      gridss: \$(java -cp "${task.ext.jarPath}" gridss.CallVariants --version 2>&1 | sed 's/-gridss//')
  END_VERSIONS
  """

  stub:
  """
  touch ${meta.id}.virusbreakend.vcf ${meta.id}.summary.tsv
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
