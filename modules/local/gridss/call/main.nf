process CALL {
  //conda (params.enable_conda ? "bioconda::gridss=2.13.2" : null)
  container 'docker.io/scwatts/gridss:2.13.2--3'

  input:
  tuple val(meta), path(bams), path(assemble_dir), val(labels)
  path gridss_config
  path genome_fa
  path genome_fai
  path genome_dict
  path genome_bwa_index_dir, stageAs: 'bwa_index'
  path genome_bwa_index_image
  path genome_gridss_index
  path blacklist

  output:
  tuple val(meta), path('gridss_call/sv_vcf.vcf.gz'), emit: vcf
  path 'versions.yml'                               , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def config_arg = gridss_config ? "--configuration ${gridss_config}" : ''
  def output_dirname = 'gridss_call'
  def labels_arg = labels.join(',')
  // NOTE(SW): Nextflow implicitly casts List<TaskPath> to an atomic TaskPath, hence the required check below
  def bams_list = bams instanceof List ? bams : [bams]
  def bams_arg = bams_list.join(' ')

  """
  # Create shadow directory with file symlinks of GRIDSS 'working' dir to prevent cache invalidation
  # NOTE: for reasons that elude me, NF doesn't always stage in the workingdir; remove if it is present
  shadow_input_directory() {
    src=\${1}
    dst="${output_dirname}/"
    for filepath_src in \$(find -L \${src} ! -type d); do
      # Get destination location for symlink
      filepath_src_rel=\$(sed 's#^'\${src}'/*##' <<< \${filepath_src})
      filepath_dst=\${dst%/}/\${filepath_src_rel}
      # Create directory for symlink
      mkdir -p \${filepath_dst%/*};
      # Get path for symlink source file, then create it
      # NOTE(SW): ideally we would get the relative path using the --relative-to but this is only
      # supported for GNU realpath and fails for others such as BusyBox, which is used in Biocontainers
      symlinkpath=\$(realpath \${filepath_src})
      ln -s \${symlinkpath} \${filepath_dst};
    done
    if [[ -L "\${src##*/}" ]]; then
      rm "\${src}"
    fi
  }
  shadow_input_directory ${assemble_dir}

  # Symlink BWA indices next to assembly FASTA
  ln -s \$(find -L ${genome_bwa_index_dir} -type f) ./

  # Run
  gridss \\
    ${args} \\
    --jvmheap ${task.memory.giga - task.ext.otherJvmHeap.giga}g \\
    --otherjvmheap ${task.ext.otherJvmHeap.giga}g \\
    --jar "${task.ext.jarPath}" \\
    --steps call \\
    --labels "${labels_arg}" \\
    --reference "${genome_fa}" \\
    --blacklist "${blacklist}" \\
    --workingdir "${output_dirname}/work/" \\
    --assembly "${output_dirname}/sv_assemblies.bam" \\
    --output "${output_dirname}/sv_vcf.vcf.gz" \\
    --threads "${task.cpus}" \\
    ${config_arg} \\
    ${bams_arg}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      gridss: \$(java -cp "${task.ext.jarPath}" gridss.CallVariants --version 2>&1 | sed 's/-gridss//')
  END_VERSIONS
  """

  stub:
  """
  mkdir -p gridss_call/
  cat <<EOF > gridss_call/sv_vcf.vcf.gz
  ##fileformat=VCFv4.1
  ##contig=<ID=.>
  #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
  .	.	.	.	.	.	.
  EOF
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
