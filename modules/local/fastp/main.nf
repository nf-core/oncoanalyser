process FASTP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--hadf994f_2' :
        'biocontainers/fastp:0.23.4--hadf994f_2' }"

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    val max_fastq_records
    val umi_location
    val umi_length
    val umi_skip

    output:
    tuple val(meta), path('*_R1.fastp.fastq.gz'), path('*_R2.fastp.fastq.gz'), emit: fastq
    path 'versions.yml'                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def split_by_lines_arg = max_fastq_records > 0 ? "--split_by_lines ${4 * max_fastq_records.toLong()}" : ''

    def umi_args_list = []
    if (umi_location) umi_args_list.add("--umi_loc ${umi_location}")
    if (umi_length) umi_args_list.add("--umi_len ${umi_length}")
    if (umi_skip >= 0) umi_args_list.add("--umi_skip ${umi_skip}")
    def umi_args = umi_args_list ? '--umi ' + umi_args_list.join(' ') : ''

    """
    fastp \\
        ${args} \\
        --in1 ${reads_fwd} \\
        --in2 ${reads_rev} \\
        ${umi_args} \\
        ${split_by_lines_arg} \\
        --thread ${task.cpus} \\
        --out1 ${meta.sample_id}_${meta.library_id}_${meta.lane}_R1.fastp.fastq.gz \\
        --out2 ${meta.sample_id}_${meta.library_id}_${meta.lane}_R2.fastp.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch 00{1..4}.${meta.sample_id}_${meta.library_id}_${meta.lane}_R1.fastp.fastq.gz
    touch 00{1..4}.${meta.sample_id}_${meta.library_id}_${meta.lane}_R2.fastp.fastq.gz

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
