process GRIDSS_PREPROCESS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sv-prep:1.2.4--hdfd78af_0' :
        'biocontainers/hmftools-sv-prep:1.2.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bam_filtered)
    path genome_fasta
    path genome_fai
    path genome_dict
    path genome_gridss_index
    path gridss_config

    output:
    tuple val(meta), path("gridss_preprocess/${meta.sample_id}.sv_prep.sorted.bam.gridss.working/"), emit: preprocess_dir
    path 'versions.yml'                                                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def config_arg = gridss_config ? "--configuration ${gridss_config}" : ''

    """
    # Symlink indices next to assembly FASTA
    ln -s \$(find -L ${genome_gridss_index} -regex '.*\\.\\(amb\\|ann\\|pac\\|gridsscache\\|sa\\|bwt\\|img\\|alt\\)') ./

    gridss_svprep \\
        ${args} \\
        --jvmheap ${Math.round(task.memory.bytes * 0.95)} \\
        --steps preprocess \\
        --reference ${genome_fasta} \\
        --workingdir gridss_preprocess/ \\
        --threads ${task.cpus} \\
        ${config_arg} \\
        --labels ${meta.sample_id} \\
        --bams ${bam} \\
        --filtered_bams ${bam_filtered} \\
        --output placeholder

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: \$(CallVariants --version 2>&1 | sed 's/-gridss\$//')
        svprep: \$(svprep -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p gridss_preprocess/${meta.sample_id}.sv_prep.sorted.bam.gridss.working/
    touch gridss_preprocess/${meta.sample_id}.sv_prep.sorted.bam.gridss.working/placeholder

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
