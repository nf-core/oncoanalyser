process FASTP {
    tag "${meta.id}"

    // TODO(MC): Resources?

    container 'docker.io/scwatts/fastp:0.23.4'

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)

    output:
    tuple val(meta), path('*_R1.fastp.fastq'), path('*_R2.fastp.fastq'), emit: fastq

    script:
    // TODO(MC): UMI flags
    //   --umi \\
    // --umi_loc per_read \\
    // --umi_len 7 \\
    // --umi_skip 1 \\

    """
    # * do not apply trimming/clipping, already done in BCL convert

    fastp \\
      --in1 ${reads_fwd} \\
      --in2 ${reads_rev} \\
      --disable_adapter_trimming \\
      --split_by_lines 40000000 \\
      --out1 ${meta.sample_id}_${meta.read_group}_R1.fastp.fastq \\
      --out2 ${meta.sample_id}_${meta.read_group}_R2.fastp.fastq
    """

    stub:
    """
    touch 00{1..4}.${meta.sample_id}_${meta.read_group}_R1.fastp.fastq
    touch 00{1..4}.${meta.sample_id}_${meta.read_group}_R2.fastp.fastq
    """
}
