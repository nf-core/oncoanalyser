process VCHORD {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker://docker.io/scwatts/hmftools-vchord:1.0.0--0'

    input:
    tuple val(meta), path(purple_dir)
    path vchord_model

    output:
    tuple val(meta), path('vchord/')
    path 'versions.yml'             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    """
    java -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        -jar /opt/vchord/vchord.jar \\
            -sample ${meta.tumor_id} \\
            -purple_dir ${purple_dir} \\
            -model ${vchord_model} \\
            -log_debug \\
            -output_dir vchord/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vchord: \$(java -jar /opt/vchord/vchord.jar -version | sed -n '/^vChord version/ { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p vchord/

    touch vchord/${meta.tumor_id}.vchord.prediction.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
