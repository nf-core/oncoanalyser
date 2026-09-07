process MINIBWA_INDEX {
    tag "$fasta"
    label 'process_single'
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    // NOTE(SW): container tag is a placeholder; finalize when container support is added
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minibwa:0.7' :
        'biocontainers/minibwa:0.7' }"

    input:
    path fasta

    output:
    path 'minibwa_index'                                   , topic: minibwa_index
    tuple val([:]), val('minibwa_index'), path('.command.*'), topic: command_files
    path 'versions.yml'                                     , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta}"
    def args = task.ext.args ?: ''

    """
    mkdir -p minibwa_index/
    minibwa index \\
        $args \\
        $fasta minibwa_index/${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minibwa: \$(minibwa version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta}"

    """
    mkdir -p minibwa_index/

    touch minibwa_index/${prefix}.l2b
    touch minibwa_index/${prefix}.mbw

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
