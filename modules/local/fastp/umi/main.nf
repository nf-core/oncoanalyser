process FASTP_UMI {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--hadf994f_2' :
        'biocontainers/fastp:0.23.4--hadf994f_2' }"

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    val umi_location
    val umi_length
    val umi_skip

    output:
    tuple val(meta), path('output/*_R1.fastp_umi.fastq.gz'), path('output/*_R2.fastp_umi.fastq.gz'), topic: fastp_umi_fastq
    tuple val(meta), val('fastp_umi'), path('.command.*')                                          , topic: command_files
    path 'versions.yml'                                                                            , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def umi_args_list = []
    if (umi_location) { umi_args_list.add("--umi_loc ${umi_location}") }
    if (umi_length) { umi_args_list.add("--umi_len ${umi_length}") }
    if (umi_skip >= 0) { umi_args_list.add("--umi_skip ${umi_skip}") }
    def umi_args = umi_args_list ? '--umi ' + umi_args_list.join(' ') : ''

    """
    mkdir -p output/

    fastp \\
        ${args} \\
        --disable_adapter_trimming \\
        --disable_length_filtering \\
        --disable_quality_filtering \\
        --disable_trim_poly_g \\
        --in1 ${reads_fwd} \\
        --in2 ${reads_rev} \\
        ${umi_args} \\
        --thread ${task.cpus} \\
        --out1 output/${meta.sample_id}_${meta.library_id}_${meta.lane}_R1.fastp_umi.fastq.gz \\
        --out2 output/${meta.sample_id}_${meta.library_id}_${meta.lane}_R2.fastp_umi.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -n '/^fastp / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    def base_name = "${meta.sample_id}_${meta.library_id}_${meta.lane}"

    """
    mkdir -p output/

    gzip <<< '' > output/${base_name}_R1.fastp_umi.fastq.gz;
    gzip <<< '' > output/${base_name}_R2.fastp_umi.fastq.gz;

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
