process NEO {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/neo:1.1_beta--0'

    input:
    val(meta)

    //tuple val(meta), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai)
    //path genome_fasta
    //val genome_ver
    //path ensembl_data_resources

    output:
    tuple val(meta), path('neo/'), emit: neo_dir
    path 'versions.yml'          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            XXX

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        neo: \$(java -jar ${task.ext.jarPath} -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p neo/
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}

