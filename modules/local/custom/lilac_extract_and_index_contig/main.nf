process EXTRACT_AND_INDEX_CONTIG {
  //conda (params.enable_conda ? "bioconda::null" : null)
  container 'docker.io/scwatts/custom-extract_and_index_contig:0.0.1--3'

  input:
  val contig_name
  path genome_fa
  path genome_fai

  output:
  path "*extracted.fa"  , emit: contig
  path "*extracted.fa.*", emit: bwa_indices
  path 'versions.yml'   , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''

  """
  samtools faidx \\
    -o ${contig_name}_extracted.fa \\
    ${genome_fa} \\
    ${contig_name}
  bwa index ${contig_name}_extracted.fa

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """

  stub:
  """
  touch ${contig_name}_extracted.fa ${contig_name}_extracted.fa.amb
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
