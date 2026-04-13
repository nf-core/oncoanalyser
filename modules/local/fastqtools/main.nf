process FASTQ_TOOLS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        '' }"

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    val max_fastq_records
    val umi_tool
    val fastp_umi_location
    val fastp_umi_length
    val fastp_umi_skip
    val fastq_tools_umi_delim
    path known_umis

    output:
    tuple val(meta),
        path('*_R1.processed.fastq.gz'),
        path('*_R2.processed.fastq.gz'), emit: fastq
    path 'versions.yml'                , emit: versions
    path '.command.*'                  , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def common_fastp_args = '--disable_quality_filtering --disable_length_filtering --disable_adapter_trimming --disable_trim_poly_g'

    //
    // UMI processing
    //
    def reads_fwd_umi = "${meta.sample_id}_${meta.library_id}_${meta.lane}_R1.umi.fastq.gz"
    def reads_rev_umi = "${meta.sample_id}_${meta.library_id}_${meta.lane}_R2.umi.fastq.gz"
    def cmd_umi_processing

    if (umi_tool == 'FASTQ_TOOLS') {

        cmd_umi_processing = """
        mkdir -p fastq_umi_processed/

        fastq-tools \\
            -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
            com.hartwig.hmftools.fastqtools.umi.FastqUmiExtracter \\
            ${args} \\
            -fastq_files '${reads_fwd};${reads_rev}' \\
            -known_umi_file ${known_umis} \\
            -umi_delim ${fastq_tools_umi_delim} \\
            -output_dir fastq_umi_processed/

        ln -sf fastq_umi_processed/*R1*fastq.gz ${reads_fwd_umi}
        ln -sf fastq_umi_processed/*R2*fastq.gz ${reads_rev_umi}
        """

    } else if (umi_tool == 'FASTP') {

        reads_fwd_umi = "${meta.sample_id}_${meta.library_id}_${meta.lane}_R1.fastq.gz"
        reads_rev_umi = "${meta.sample_id}_${meta.library_id}_${meta.lane}_R2.fastq.gz"

        def umi_args_list = []
        umi_args_list.add('--umi')
        if (umi_location)  { umi_args_list.add("--umi_loc ${umi_location}") }
        if (umi_length)    { umi_args_list.add("--umi_len ${umi_length}") }
        if (umi_skip >= 0) { umi_args_list.add("--umi_skip ${umi_skip}") }
        def umi_args = umi_args_list.join(' ')

        cmd_umi_processing = """
        fastp \\
            ${common_fastp_args} \\
            --in1 ${reads_fwd} \\
            --in2 ${reads_rev} \\
            ${umi_args} \\
            --thread ${task.cpus} \\
            --out1 ${reads_fwd_umi} \\
            --out2 ${reads_rev_umi}
        """
        
    } else {
    
        cmd_umi_processing = ''
        reads_fwd_umi = reads_fwd
        reads_rev_umi = reads_rev
        
    }

    //
    // FASTQ splitting
    //
    def reads_fwd_split = "${meta.sample_id}_${meta.library_id}_${meta.lane}_R1.processed.fastq.gz"
    def reads_rev_split = "${meta.sample_id}_${meta.library_id}_${meta.lane}_R2.processed.fastq.gz"
    def cmd_fastq_split
    
    if (max_fastq_records > 0) {

        cmd_fastq_split = """
        fastp \\
            ${common_fastp_args} \\
            --in1 ${reads_fwd_umi} \\
            --in2 ${reads_rev_umi} \\
            --split_by_lines ${4 * max_fastq_records.toLong()} \\
            --thread ${task.cpus} \\
            --out1 ${reads_fwd_split} \\
            --out2 ${reads_rev_split}
        """

    } else {

        cmd_fastq_split = """
        ln -sf ${reads_fwd_umi} ${reads_fwd_split}
        ln -sf ${reads_rev_umi} ${reads_rev_split}
        """

    }

    """
    ${cmd_umi_processing}
    
    ${cmd_fastq_split}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/^.* //')
        fastq-tools: 1.0
    END_VERSIONS
    """

    stub:
    """
    touch output/00{1..4}.${meta.sample_id}_${meta.library_id}_${meta.lane}_R1.fastq.gz
    touch output/00{1..4}.${meta.sample_id}_${meta.library_id}_${meta.lane}_R2.fastq.gz

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
