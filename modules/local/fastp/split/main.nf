nextflow.enable.types = true

process FASTP_SPLIT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--hadf994f_2' :
        'biocontainers/fastp:0.23.4--hadf994f_2' }"

    input:
    tuple(meta: Record, reads_fwd: Path, reads_rev: Path)
    max_fastq_records: Integer

    topic:
    tuple(meta, files('output/*_R1.fastp_split.fastq.gz'), files('output/*_R2.fastp_split.fastq.gz', optional: true)) >> 'fastp_split_fastq'
    tuple(meta, 'fastp_split', files('.command.*'))                                                                   >> 'command_files'
    file('versions.yml')                                                                                              >> 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def in2_args = meta.single_end ? '' : "--in2 ${reads_rev}"
    def out2_args = meta.single_end ? '' : "--out2 output/${meta.sample_id}_${meta.library_id}_${meta.lane}_R2.fastp_split.fastq.gz"

    """
    mkdir -p output/

    fastp \\
        ${args} \\
        --disable_adapter_trimming \\
        --disable_length_filtering \\
        --disable_quality_filtering \\
        --disable_trim_poly_g \\
        --in1 ${reads_fwd} \\
        ${in2_args} \\
        --split_by_lines ${4 * max_fastq_records.toLong()} \\
        --thread ${task.cpus} \\
        --out1 output/${meta.sample_id}_${meta.library_id}_${meta.lane}_R1.fastp_split.fastq.gz \\
        ${out2_args}

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
