process FASTP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--hadf994f_2' :
        'biocontainers/fastp:0.23.4--hadf994f_2' }"

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    val(max_fastq_records)
    val umi_length

    output:
    tuple val(meta), path('*_R1.fastp.fastq.gz'), path('*_R2.fastp.fastq.gz'), emit: fastq
    path 'versions.yml'                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def umi_extraction = umi_length > 0 ? '--umi --umi_loc per_read --umi_len ' + umi_length : ''

    """
    # * do not apply trimming/clipping, already done in BCL convert
    # * turn off all filtering
    # * extract UMIs if a UMI length has been specified

    """
    fastp \\
        ${args} \\
        --in1 ${reads_fwd} \\
        --in2 ${reads_rev} \\
        --disable_quality_filtering \\
        --disable_length_filtering \\
        --disable_adapter_trimming \\
        --disable_trim_poly_g \\
        --split_by_lines ${4 * max_fastq_records} \\
        ${umi_extraction} \\
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
