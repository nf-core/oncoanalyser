process VIRUSINTERPRETER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-virus-interpreter:1.7--hdfd78af_0' :
        'biocontainers/hmftools-virus-interpreter:1.7--hdfd78af_0' }"

    input:
    tuple val(meta), path(virus_tsv), path(purple_dir), path(bamtools_somatic_dir)
    path taxonomy_db
    path reporting_db
    path blocklist_db

    output:
    tuple val(meta), path('virusinterpreter/'), emit: virusinterpreter_dir
    path 'versions.yml'                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    """
    mkdir -p virusinterpreter/

    virusinterpreter \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -purple_dir ${purple_dir} \\
        -tumor_metrics_dir ${bamtools_somatic_dir} \\
        -virus_breakend_tsv ${virus_tsv} \\
        -taxonomy_db_tsv ${taxonomy_db} \\
        -virus_reporting_db_tsv ${reporting_db} \\
        -virus_blacklisting_db_tsv ${blocklist_db} \\
        -output_dir virusinterpreter/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        virusinterpreter: \$(virusinterpreter -version | sed -n '/^VirusInterpreter version/ { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p virusinterpreter/
    touch virusinterpreter/${meta.sample_id}.virus.annotated.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
