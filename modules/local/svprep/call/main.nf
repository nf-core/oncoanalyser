process GRIDSS_CALL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sv-prep:1.2.3--hdfd78af_1' :
        'quay.io/biocontainers/hmftools-sv-prep:1.2.3--hdfd78af_1' }"

    input:
    tuple val(meta), path(bams), path(bams_filtered), path(assemble_dir), val(labels)
    path gridss_config
    path genome_fasta
    path genome_fai
    path genome_dict
    path genome_bwa_index_dir, stageAs: 'bwa_index'
    path genome_bwa_index_image
    path genome_gridss_index
    path blocklist

    output:
    tuple val(meta), path('gridss_call/sv.svprep.gridss.vcf.gz'), emit: vcf
    path 'versions.yml'                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config_arg = gridss_config ? "--configuration ${gridss_config}" : ''
    def output_dirname = 'gridss_call'
    def labels_list = labels instanceof List ? labels : [labels]
    def labels_arg = labels_list.join(',')
    // NOTE(SW): Nextflow implicitly casts List<TaskPath> to an atomic TaskPath, hence the required check below
    def bams_list = bams instanceof List ? bams : [bams]
    def bams_arg = "--bams ${bams_list.join(',')}"
    def bams_filtered_list = bams_filtered instanceof List ? bams_filtered : [bams_filtered]
    def bams_filtered_arg = "--filtered_bams ${bams_filtered_list.join(',')}"

    // JVM heap for other tasks must be no greater than 1/4 of task memory, defaults to 1 GB if not provided
    def otherJvmHeap = Math.min(
        Math.round(task.memory.bytes * 0.25),
        task.ext.otherJvmHeap ? task.ext.otherJvmHeap.bytes : 1.GB.bytes
    )

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
    gridss_svprep \\
        ${args} \\
        --jvmheap ${Math.round((task.memory.bytes - otherJvmHeap) * 0.95)} \\
        --otherjvmheap ${otherJvmHeap} \\
        --steps call \\
        --labels ${labels_arg} \\
        --reference ${genome_fasta} \\
        --blacklist ${blocklist} \\
        --workingdir ${output_dirname}/work/ \\
        --assembly ${output_dirname}/sv_assemblies.bam \\
        --output ${output_dirname}/sv.svprep.gridss.vcf.gz \\
        --threads ${task.cpus} \\
        ${config_arg} \\
        ${bams_arg} \\
        ${bams_filtered_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: \$(CallVariants --version 2>&1 | sed 's/-gridss\$//')
        svprep: \$(svprep -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p gridss_call/
    cat <<EOF > gridss_call/sv.svprep.gridss.vcf.gz
    ##fileformat=VCFv4.1
    ##contig=<ID=.>
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    .	.	.	.	.	.	.
    EOF
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
