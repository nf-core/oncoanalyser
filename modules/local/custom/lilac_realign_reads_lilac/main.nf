process REALIGN_READS {
  //conda (params.enable_conda ? "bioconda::hmftools-lilac=1.2" : null)
  container 'docker.io/scwatts/custom-realign_reads_lilac:0.0.1--3'

  input:
  tuple val(meta), path(bam), path(bai)
  path reference
  path reference_indices

  output:
  tuple val(meta), path("*realigned.bam"), path("*realigned.bam.bai"), emit: bam
  path 'versions.yml'                                                , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''

  """
  sambamba sort -n ${bam} -o ${meta.id}_sorted.bam
  samtools fastq -@${task.threads} ${meta.id}_sorted.bam \\
      -1 ${meta.id}_R1.fastq.gz \\
      -2 ${meta.id}_R2.fastq.gz \\
      -0 ${meta.id}_other.fastq.gz \\
      -s ${meta.id}_singleton.fastq.gz;

  bwa mem \\
    -t${task.cpus} \\
    -Y \\
    ${reference} \\
    ${meta.id}_R1.fastq.gz \\
    ${meta.id}_R2.fastq.gz | \\
    samtools sort -T tmp -o ${meta.id}_realigned.bam
  samtools index ${meta.id}_realigned.bam

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
      sambamba: \$(sambamba --version 2>&1 | sed -n '/sambamba/ s/^sambamba \\(.\\+\\)/\\1/p' | head -n1)
  END_VERSIONS
  """

  stub:
  """
  touch ${meta.id}_realigned.bam ${meta.id}_realigned.bam.bai
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
