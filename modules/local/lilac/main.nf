process LILAC {
  //conda (params.enable_conda ? "bioconda::hmftools-lilac=1.2" : null)
  container 'docker.io/scwatts/lilac:1.4.1--0'

  input:
  tuple val(meta), path(tumor_bam), path(normal_bam), path(tumour_bai), path(normal_bai), path(purple_dir)
  path genome_fa
  val genome_ver
  path lilac_resource_dir, stageAs: 'lilac_resource_dir'

  output:
  tuple val(meta), path('lilac/'), emit: lilac_dir
  path 'versions.yml'            , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def sample_name = get_sample_name(meta, tumor_bam, normal_bam)
  def tumor_bam_arg = tumor_bam ? "-tumor_bam ${tumor_bam}" : ''
  def purple_args = purple_dir ? """
    -gene_copy_number ${purple_dir}/${sample_name}.purple.cnv.gene.tsv \\
    -somatic_vcf ${purple_dir}/${sample_name}.purple.sv.vcf.gz \\
  """ : ''

  """
  java \\
    -Xmx${task.memory.giga}g \\
    -jar "${task.ext.jarPath}" \\
      ${args} \\
      -sample "${sample_name}" \\
      ${tumor_bam_arg} \\
      -reference_bam ${normal_bam} \\
      -ref_genome_version "${genome_ver}" \\
      -ref_genome "${genome_fa}" \\
      -resource_dir "${lilac_resource_dir}" \\
      ${purple_args.replaceAll('\\n', '')} \\
      -threads "${task.cpus}" \\
      -output_dir lilac/

  # NOTE(SW): hard coded since there is no reliable way to obtain version information.
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      lilac: 1.4.1
  END_VERSIONS
  """

  stub:
  """
  mkdir -p lilac/
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}

def get_sample_name(meta, tumor_bam, normal_bam) {
    if (tumor_bam) {
        return meta.get(['sample_name', 'tumor'])
    } else if (normal_bam) {
        return meta.get(['sample_name', 'normal'])
    } else {
        Sys.exit(1)
    }
}
