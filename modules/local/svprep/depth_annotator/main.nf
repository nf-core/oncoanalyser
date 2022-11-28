process DEPTH_ANNOTATOR {
  container 'docker.io/scwatts/svprep:1.1--0'

  input:
  tuple val(meta), path(bams), path(bais), path(vcf), val(labels)
  path genome_fa
  val genome_ver

  output:
  tuple val(meta), path("*.vcf.gz"), emit: vcf
  path 'versions.yml'              , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def labels_arg = labels.join(',')
  // NOTE(SW): Nextflow implicitly casts List<TaskPath> to an atomic TaskPath, hence the required check below
  def bams_list = bams instanceof List ? bams : [bams]
  def bams_arg = "${bams_list.join(',')}"

  """
  java \\
    -Xmx${task.memory.giga}g \\
    -cp "${task.ext.jarPath}" com.hartwig.hmftools.svprep.depth.DepthAnnotator \\
      ${args} \\
      -input_vcf ${vcf} \\
      -output_vcf sv.svprep.gridss.depths.vcf.gz \\
      -samples "${labels_arg}" \\
      -bam_files "${bams_arg}" \\
      -ref_genome ${genome_fa} \\
      -ref_genome_version ${genome_ver} \\
      -threads ${task.cpus}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      svprep: \$(java -jar "${task.ext.jarPath}" 2>&1 | head -n1 | sed 's/^.*SvPrep version: //')
  END_VERSIONS
  """

  stub:
  """
  touch sv_vcf.depths.vcf.gz
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
