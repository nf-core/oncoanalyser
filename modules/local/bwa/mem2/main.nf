process BWA_MEM2 {
    tag "${meta.id}"
    label 'process_high'

    // TODO(SW): create container
    //container 'foo'

    input:
    // TODO(SW): decide input structure
    tuple val(meta), path(fastqs)

    output:
    // TODO(SW): set outputs
    tuple val(meta), path('bar'), emit: bam
    path 'versions.yml'         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    // TODO(SW): implement process
    """
    touch bar

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: foo
    END_VERSIONS
    """

    stub:
    // TODO(SW): implement stub
    """
    touch bar
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
