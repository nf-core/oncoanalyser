process COBALT {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/cobalt:1.13--1'

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai)
    path gc_profile

    output:
    tuple val(meta), path('cobalt/'), emit: cobalt_dir
    path 'versions.yml'             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    java \\
        -Xmx${task.memory.giga}g \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -tumor ${meta.tumor_id} \\
            -tumor_bam ${tumor_bam} \\
            -reference ${meta.normal_id} \\
            -reference_bam ${normal_bam} \\
            -output_dir cobalt/ \\
            -threads ${task.cpus} \\
            -gc_profile ${gc_profile}

    # NOTE(SW): hard coded since there is no reliable way to obtain version information.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobalt: 1.13
    END_VERSIONS
    """

    stub:
    """
    mkdir -p cobalt/
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
