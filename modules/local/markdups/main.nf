process MARKDUPS {
    tag "${meta.id}"
    label 'process_low'

    // TODO(SW): create container
    //container 'foo'

    input:
    tuple val(meta), path(bam)

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
    echo bar

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        markdups: foo
    END_VERSIONS
    """

    stub:
    // TODO(SW): implement stub
    """
    touch bar
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
