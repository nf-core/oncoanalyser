process VIRUSINTERPRETER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-virus-interpreter:1.7.2--hdfd78af_0' :
        'biocontainers/hmftools-virus-interpreter:1.7.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(virusbreakend_tsv), path(bamtools_somatic_dir), path(purple_dir)
    path taxonomy_db
    path reporting_db
    path blocklist_db

    output:
    tuple val(meta), path('virusinterpreter/')                  , topic: virusinterpreter_dir
    tuple val(meta), val('virusinterpreter'), path('.command.*'), topic: command_files
    path 'versions.yml'                                         , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    """
    mkdir -p virusinterpreter/

    virusinterpreter \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -purple_dir ${purple_dir} \\
        -tumor_metrics_dir ${bamtools_somatic_dir} \\
        -virus_breakend_tsv ${virusbreakend_tsv} \\
        -taxonomy_db_tsv ${taxonomy_db} \\
        -virus_reporting_db_tsv ${reporting_db} \\
        -virus_blacklisting_db_tsv ${blocklist_db} \\
        ${log_level_arg} \\
        -output_dir virusinterpreter/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        virusinterpreter: \$(virusinterpreter -version | sed -n '/^VirusInterpreter version/ { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p virusinterpreter/

    touch virusinterpreter/${meta.sample_id}.virus.annotated.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
