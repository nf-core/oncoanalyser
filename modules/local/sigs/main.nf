process SIGS {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/sigs:1.0--0'

    input:
    tuple val(meta), path(smlv_vcf)
    path signatures

    output:
    tuple val(meta), path('sigs/'), emit: sigs_dir
    path 'versions.yml'           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p sigs/

    java \\
        -Xmx${task.memory.giga}g \\
        -jar ${task.ext.jarPath} \\
            -sample ${meta.id} \\
            -somatic_vcf_file ${smlv_vcf} \\
            -signatures_file ${signatures} \\
            -output_dir ./sigs/

    # NOTE(SW): hard coded since there is no reliable way to obtain version information.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sigs: 1.0
    END_VERSIONS
    """

    stub:
    """
    mkdir sigs/
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
