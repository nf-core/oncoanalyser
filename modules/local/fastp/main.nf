process FASTP {
    tag "${meta.subject_id}__${meta.sample_id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--hadf994f_2' :
        'quay.io/biocontainers/fastp:0.23.4--hadf994f_2' }"

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    val(max_fastq_records)

    output:
    tuple val(meta), path('*_R1.fastp.fastq'), path('*_R2.fastp.fastq'), emit: fastq
    path 'versions.yml'                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # * do not apply trimming/clipping, already done in BCL convert
    # * turn off all filtering
    # * do not process umis, already done for us

    fastp \\
      --in1 ${reads_fwd} \\
      --in2 ${reads_rev} \\
      --disable_quality_filtering \\
      --disable_length_filtering \\
      --disable_adapter_trimming \\
      --disable_trim_poly_g \\
      --split_by_lines ${4 * max_fastq_records} \\
      --out1 ${meta.sample_id}_${meta.read_group}_R1.fastp.fastq \\
      --out2 ${meta.sample_id}_${meta.read_group}_R2.fastp.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch 00{1..4}.${meta.sample_id}_${meta.read_group}_R1.fastp.fastq
    touch 00{1..4}.${meta.sample_id}_${meta.read_group}_R2.fastp.fastq

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
