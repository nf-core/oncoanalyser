process TAUR {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-taur:1.0--hdfd78af_0' :
        'biocontainers/hmftools-taur:1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    val umi_delim
    path known_umis

    output:
    tuple val(meta), path('output/*R1.umi.fastq.gz'), path('output/*R2.umi.fastq.gz'), topic: taur_fastq
    tuple val(meta), val('taur'), path('.command.*')                                 , topic: command_files
    path 'versions.yml'                                                              , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    """
    mkdir -p output/

    taur \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        com.hartwig.hmftools.taur.TaurApplication \\
        ${args} \\
        -fastq_files '${reads_fwd};${reads_rev}' \\
        -known_umi_file ${known_umis} \\
        -umi_delim ${umi_delim} \\
        ${log_level_arg} \\
        -output_dir output/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taur: 1.0
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p output/

    gzip <<< '' > output/${meta.sample_id}_${meta.library_id}_${meta.lane}_R1.umi.fastq.gz
    gzip <<< '' > output/${meta.sample_id}_${meta.library_id}_${meta.lane}_R2.umi.fastq.gz

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
