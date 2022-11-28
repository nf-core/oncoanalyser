process SVPREP {
  container 'docker.io/scwatts/svprep:1.1--0'

  input:
  tuple val(meta), path(bam), path(bai), path(junctions)
  path genome_fa
  val genome_ver
  path sv_blacklist
  path known_fusions
  val write_types
  val calc_fragment_length

  output:
  tuple val(meta), path("*.sorted.bam")           , emit: bam
  tuple val(meta), path("*.sv_prep.junctions.csv"), emit: junctions
  path 'versions.yml'                             , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def write_types_arg = write_types ? "-write_types \'${write_types}\'" : ''
  def existing_juction_file_arg = junctions ? "-existing_junction_file ${junctions}" : ""
  def calc_fragment_length_arg = calc_fragment_length ? "-calc_fragment_length" : ""

  """
  java \\
    -Xmx${task.memory.giga}g \\
    -jar "${task.ext.jarPath}" \\
      ${args} \\
      -sample "${meta.id}" \\
      -bam_file "${bam}" \\
      -ref_genome "${genome_fa}" \\
      -ref_genome_version "${genome_ver}" \\
      -blacklist_bed "${sv_blacklist}" \\
      -known_fusion_bed "${known_fusions}" \\
      ${write_types_arg} \\
      ${existing_juction_file_arg} \\
      ${calc_fragment_length_arg} \\
      -threads "${task.cpus}" \\
      -output_dir ./

  samtools sort -O bam "${meta.id}.sv_prep.bam" -o "${meta.id}.sv_prep.sorted.bam"

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      sage: \$(java -jar "${task.ext.jarPath}" 2>&1 | head -n1 | sed 's/^.*SvPrep version: //')
  END_VERSIONS
  """

  stub:
  """
  touch "${meta.id}.sv_prep.sorted.bam"
  touch "${meta.id}.sv_prep.junctions.csv"
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
