process CHORD {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/chord:2.00--0'

    input:
    tuple val(meta), path(smlv_vcf), path(sv_vcf)
    val genome_ver

    output:
    path '*_chord_signatures.txt', emit: signatures
    path '*_chord_prediction.txt', emit: prediction
    path 'versions.yml'          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    extractSigPredictHRD.R \\
        ./ \\
        ${meta.id} \\
        ${smlv_vcf} \\
        ${sv_vcf} \\
        ${genome_ver} \\
        ${meta.id}_chord_signatures.txt \\
        ${meta.id}_chord_prediction.txt

    # NOTE(SW): hard coded since there is no reliable way to obtain version information.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CHORD: \$(R -s -e "message(packageVersion('CHORD'))")
        mutSigExtractor: \$(R -s -e "message(packageVersion('mutSigExtractor'))")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_chord_signatures.txt
    touch ${meta.id}_chord_prediction.txt
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
