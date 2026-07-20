process SEQKIT_SPLIT2 {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.13.0--he881be0_0' :
        'biocontainers/seqkit:2.13.0--he881be0_0' }"

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    val max_fastq_records

    output:
    tuple val(meta), path('output/*_R1.*.fastq.gz'), path('output/*_R2*.*.fastq.gz'), topic: seqkit_split2_fastq
    tuple val(meta), val('seqkit_split2'), path('.command.*')                       , topic: command_files
    path 'versions.yml'                                                             , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    seqkit split2 \\
        --by-size ${max_fastq_records} \\
        --read1 ${reads_fwd} \\
        --read2 ${reads_rev} \\
        --threads ${task.cpus} \\
        --by-size-prefix ${meta.sample_id}_${meta.library_id}_${meta.lane}_R{read}.split_ \\
        --extension .gz \\
        --out-dir output/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed -n '/^seqkit/ { s/^seqkit v//p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p output/

    touch output/${meta.sample_id}_${meta.library_id}_${meta.lane}_R1.split_00{1..4}.fastq.gz
    touch output/${meta.sample_id}_${meta.library_id}_${meta.lane}_R2.split_00{1..4}.fastq.gz

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
