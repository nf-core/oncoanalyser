process PREPROCESS {
  container 'docker.io/scwatts/svprep:1.1--0'

  input:
  tuple val(meta), path(bam), path(bam_filtered)
  path gridss_config
  path genome_fa
  path genome_fai
  path genome_dict
  path genome_bwa_index_dir, stageAs: 'bwa_index'
  path genome_bwa_index_image
  path genome_gridss_index

  output:
  tuple val(meta), path("gridss_preprocess/${meta.id}.sv_prep.sorted.bam.gridss.working/"), emit: preprocess_dir
  path 'versions.yml'                                                                     , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def config_arg = gridss_config ? "--configuration ${gridss_config}" : ''

  """
  # Symlink BWA indices next to assembly FASTA
  ln -s \$(find -L ${genome_bwa_index_dir} -type f) ./

  gridss_svprep \\
    ${args} \\
    --jvmheap "${task.memory.giga}g" \\
    --jar "${task.ext.jarPath}" \\
    --steps preprocess \\
    --reference "${genome_fa}" \\
    --workingdir gridss_preprocess/ \\
    --threads "${task.cpus}" \\
    ${config_arg} \\
    --labels ${meta.id} \\
    --bams ${bam} \\
    --filtered_bams ${bam_filtered} \\
    --output placeholder

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      gridss: \$(java -cp "${task.ext.jarPath}" gridss.CallVariants --version 2>&1 | sed 's/-gridss//')
  END_VERSIONS
  """

  stub:
  """
  mkdir -p gridss_preprocess/${meta.id}.sv_prep.sorted.bam.gridss.working/
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
