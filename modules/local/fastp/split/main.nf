process FASTP_SPLIT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--hadf994f_2' :
        'biocontainers/fastp:0.23.4--hadf994f_2' }"

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    val max_fastq_records

    output:
    tuple val(meta), path('output/*_R1.fastp_split.fastq.gz'), path('output/*_R2.fastp_split.fastq.gz'), topic: fastp_split_fastq
    tuple val(meta), val('fastp_split'), path('.command.*')                                            , topic: command_files
    path 'versions.yml'                                                                                , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

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
        --split_by_lines ${4 * max_fastq_records.toLong()} \\
        --thread ${task.cpus} \\
        --out1 output/${meta.sample_id}_${meta.library_id}_${meta.lane}_R1.fastp_split.fastq.gz \\
        --out2 output/${meta.sample_id}_${meta.library_id}_${meta.lane}_R2.fastp_split.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    def base_name = "${meta.sample_id}_${meta.library_id}_${meta.lane}"

    """
    mkdir -p output/

    gzip <<< '' > output/001.${base_name}_R1.fastp_split.fastq.gz;
    gzip <<< '' > output/002.${base_name}_R1.fastp_split.fastq.gz;
    gzip <<< '' > output/003.${base_name}_R1.fastp_split.fastq.gz;
    gzip <<< '' > output/004.${base_name}_R1.fastp_split.fastq.gz;

    gzip <<< '' > output/001.${base_name}_R2.fastp_split.fastq.gz;
    gzip <<< '' > output/002.${base_name}_R2.fastp_split.fastq.gz;
    gzip <<< '' > output/003.${base_name}_R2.fastp_split.fastq.gz;
    gzip <<< '' > output/004.${base_name}_R2.fastp_split.fastq.gz;

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
