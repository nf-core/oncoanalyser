process VCUPPA {
    tag "${meta.id}"
    label 'process_low'

    container 'docker://docker.io/scwatts/hmftools-vcuppa:0.1.0--0'

    input:
    tuple val(meta), path(purple_dir)
    val genome_ver
    path vcuppa_model
    path vcuppa_features_list

    output:
    tuple val(meta), path('vcuppa/')
    path 'versions.yml'             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    """
    # Extract input features
    mkdir -p vcuppa/

    cuppa \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        com.hartwig.hmftools.cup.prep.CuppaDataPrep \\
        ${args} \\
        -sample ${meta.tumor_id} \\
        -categories DNA \\
        -purple_dir ${purple_dir} \\
        -ref_genome_version ${genome_ver} \\
        -output_dir vcuppa/

    # Run prediction
    java -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        -jar /opt/vcuppa/vcuppa.jar \\
          ${args2} \\
          -sample ${meta.tumor_id} \\
          -purple_dir ${purple_dir} \\
          -cuppa_data vcuppa/${meta.tumor_id}.cuppa_data.tsv.gz \\
          -model ${vcuppa_model} \\
          -cuppa_data_definitions ${vcuppa_features_list} \\
          -log_debug \\
          -output_dir vcuppa/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cuppa: \$(cuppa -version | sed -n '/Cuppa version/ { s/^.* //p }')
        vcuppa: \$(java -jar /opt/vcuppa/vcuppa.jar -version | sed -n '/vCuppa version/ { s/^.* //p }'
    END_VERSIONS
    """

    stub:
    """
    mkdir -p vcuppa/

    touch vcuppa/${meta.tumor_id}.cuppa_data.tsv.gz
    touch vcuppa/${meta.tumor_id}.vcuppa.prediction.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
