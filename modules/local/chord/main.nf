process CHORD {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/chord:2.00--0'

    input:
    tuple val(meta), path(smlv_vcf), path(sv_vcf)
    val genome_ver

    output:
    tuple val(meta), path('chord/'), emit: chord_dir
    path 'versions.yml'            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p chord/

    extractSigPredictHRD.R \\
        ./ \\
        ${meta.sample_id} \\
        ${smlv_vcf} \\
        ${sv_vcf} \\
        ${genome_ver} \\
        chord_signatures.txt \\
        chord_prediction.txt

    mv ${meta.id}_chord_signatures.txt ${meta.id}_chord_prediction.txt chord/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CHORD: \$(R -s -e "message(packageVersion('CHORD'))")
        mutSigExtractor: \$(R -s -e "message(packageVersion('mutSigExtractor'))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p chord/
    touch chord/${meta.sample_id}_chord_signatures.txt
    touch chord/${meta.sample_id}_chord_prediction.txt
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
