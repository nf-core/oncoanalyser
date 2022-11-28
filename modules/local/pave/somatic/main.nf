// NOTE(SW): PAVE gnomad filtering is not yet documented but is used in Pipeline5 https://github.com/hartwigmedical/pipeline5/blob/master/cluster/src/main/java/com/hartwig/pipeline/tertiary/pave/PaveArguments.java#L27-L28

process PAVE_SOMATIC {
  //conda (params.enable_conda ? "bioconda::hmftools-pave=1.3" : null)
  container 'docker.io/scwatts/pave:1.4--0'

  input:
  tuple val(meta), path(sage_vcf)
  path genome_fa
  path genome_fai
  val genome_ver
  path sage_pon_file
  path mappability_bed
  path driver_gene_panel
  path ensembl_data_dir

  output:
  tuple val(meta), path("*.vcf.gz")    , emit: vcf
  tuple val(meta), path("*.vcf.gz.tbi"), emit: index
  path 'versions.yml'                  , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def pon_filters = "HOTSPOT:5:5;PANEL:2:5;UNKNOWN:2:0"

  """
  java \\
    -Xmx${task.memory.giga}g \\
    -jar "${task.ext.jarPath}" \\
      -sample "${meta.get(['sample_name', 'tumor'])}" \\
      -ref_genome_version "${genome_ver}" \\
      -ref_genome "${genome_fa}" \\
      -ensembl_data_dir "${ensembl_data_dir}" \\
      -driver_gene_panel "${driver_gene_panel}" \\
      -pon_file "${sage_pon_file}" \\
      -pon_filters "${pon_filters}" \\
      -mappability_bed "${mappability_bed}" \\
      -vcf_file "${sage_vcf}" \\
      -read_pass_only \\
      -output_dir ./

  # NOTE(SW): hard coded since there is no reliable way to obtain version information.
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      pave: 1.4
  END_VERSIONS
  """

  stub:
  """
  touch ${meta.get(['sample_name', 'tumor'])}.sage.pave_somatic.vcf.gz{,.tbi}
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
