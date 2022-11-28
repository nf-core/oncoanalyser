process SAGE_SOMATIC {
  //conda (params.enable_conda ? "bioconda::hmftools-sage=3.1" : null)
  container 'docker.io/scwatts/sage:3.2--0'

  input:
  tuple val(meta), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai)
  path genome_fa
  path genome_fai
  path genome_dict
  val genome_ver
  path sage_known_hotspots_somatic
  path sage_coding_panel
  path sage_high_confidence
  path ensembl_data_dir

  output:
  tuple val(meta), path("${meta.subject_name}.sage_somatic.vcf.gz"), emit: vcf
  path 'versions.yml'                                              , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''

  """
  java \\
    -Xmx${task.memory.giga}g \\
    -jar "${task.ext.jarPath}" \\
      ${args} \\
      -reference "${meta.get(['sample_name', 'normal'])}" \\
      -reference_bam "${normal_bam}" \\
      -tumor "${meta.get(['sample_name', 'tumor'])}" \\
      -tumor_bam "${tumor_bam}" \\
      -ref_genome_version "${genome_ver}" \\
      -ref_genome "${genome_fa}" \\
      -hotspots "${sage_known_hotspots_somatic}" \\
      -panel_bed "${sage_coding_panel}" \\
      -high_confidence_bed "${sage_high_confidence}" \\
      -ensembl_data_dir "${ensembl_data_dir}" \\
      -threads "${task.cpus}" \\
      -out "${meta.subject_name}.sage_somatic.vcf.gz"

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      sage: \$(java -jar "${task.ext.jarPath}" | head -n1 | sed 's/.*Sage version: //')
  END_VERSIONS
  """

  stub:
  """
  touch "${meta.subject_name}.sage_somatic.vcf.gz"
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
