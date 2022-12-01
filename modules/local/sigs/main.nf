process SIGS {
    container 'docker.io/scwatts/sigs:1.1--0'

    input:
    tuple val(meta), path(smlv_vcf)
    path signatures_file

    output:
    tuple val(meta), path('sigs/'), emit: sigs_dir
    path 'versions.yml'           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    java \\
        -Xmx${task.memory.giga}g \\
        -jar ${task.ext.jarPath} \\
            -sample ${meta.id} \\
            -signatures_file ${signatures_file} \\
            -output_dir ./sigs/

    # NOTE(SW): hard coded since there is no reliable way to obtain version information.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sigs: 1.1
    END_VERSIONS
    """

    stub:
    """
    mkdir sigs/
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
