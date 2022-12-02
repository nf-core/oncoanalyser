process CUPPA_VISUALISER {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/cuppa:1.7--0'

    input:
    tuple val(meta), path(cuppa_csv)

    output:
    path '*png'
    path '*txt'
    path 'versions.yml'    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    python ${task.ext.pythonPath} \\
        -sample ${meta.id} \\
        -sample_data ${cuppa_csv} \\
        -output_dir ./

    # NOTE(SW): hard coded since there is no reliable way to obtain version information.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cuppa: 1.7
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.cuppa.chart.png
    touch ${meta.id}.cuppa.conclusion.txt
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
