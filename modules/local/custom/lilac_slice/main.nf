process SLICE {
  //conda (params.enable_conda ? "bioconda::samtools:1.15.1" : null)
  container 'quay.io/biocontainers/samtools:1.15.1--h1170115_0'

  input:
  tuple val(meta), path(bam), path(bai), path(bed)

  output:
  tuple val(meta), path("*sliced.bam"), path("*sliced.bam.bai"), emit: bam
  path 'versions.yml'                                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''

  """
  samtools view \\
    ${args} \\
    -L "${bed}" \\
    -@${task.cpus} \\
    -Obam \\
    ${bam} | \\
    samtools sort -T tmp -o "${bam.simpleName}.sliced.bam"

  samtools index "${bam.simpleName}.sliced.bam"

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """

  stub:
  """
  touch ${bam.simpleName}.sliced.bam ${bam.simpleName}.sliced.bam.bai
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
