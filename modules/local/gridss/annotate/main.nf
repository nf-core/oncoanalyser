process GRIDSS_ANNOTATE {
  tag "${meta.id}"
  label 'process_single'

  container 'docker.io/scwatts/gridss:2.13.2--3'

  input:
  tuple val(meta), path(gridss_vcf)

  output:
  tuple val(meta), path('gridss_annotate/*.gridss.annotated.vcf.gz'), emit: vcf
  path 'versions.yml'                                               , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''

  """
  gridss_annotate_vcf_repeatmasker \\
    ${args} \\
    --jar ${task.ext.jarPath} \\
    --output gridss_annotate/${meta.id}.gridss.annotated.vcf.gz \\
    --workingdir gridss_annotate/work/ \\
    --threads ${task.cpus} \\
    ${gridss_vcf}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      gridss: \$(java -cp ${task.ext.jarPath} gridss.CallVariants --version 2>&1 | sed 's/-gridss//')
  END_VERSIONS
  """

  stub:
  """
  mkdir -p gridss_annotate/
  touch gridss_annotate/${meta.id}.gridss.annotated.vcf.gz
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
