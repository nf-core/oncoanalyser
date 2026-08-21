nextflow.enable.types = true

process VIRUSINTERPRETER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-virus-interpreter:1.7.2--hdfd78af_0' :
        'biocontainers/hmftools-virus-interpreter:1.7.2--hdfd78af_0' }"

    input:
    tuple(meta: Map, virusbreakend_tsv: Path, bamtools_somatic_dir: Path, purple_dir: Path)
    taxonomy_db: Path
    reporting_db: Path
    blocklist_db: Path

    topic:
    tuple(meta, file('virusinterpreter/'))               >> 'virusinterpreter_dir'
    tuple(meta, 'virusinterpreter', files('.command.*')) >> 'command_files'
    file('versions.yml')                                 >> 'versions'

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
