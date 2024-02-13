process BWA_MEM2 {
    tag "${meta.id}"
    label 'process_high'

    // TODO(MC): Upload container.
    container 'bwa-mem2:2.2.1-sambamba'

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    path genome_fasta
    // TODO(MC): Copied into local genome_bwa_index for ref genome 37:
    //    + Homo_sapiens.GRCh37.GATK.illumina.fasta.bwt.2bit.64
    //    + Homo_sapiens.GRCh37.GATK.illumina.fasta.0123
    path genome_bwa_index

    output:
    tuple val(meta), path('*.bam'), emit: bam
    // TODO(MC): Versions.
    // path 'versions.yml'         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // # TODO(MC): read group
    // # -R ${meta.read_group}

    """
    ln -s \$(find -L ${genome_bwa_index} -type f) ./

    bwa-mem2 mem \\
    -Y \\
    -t ${task.cpus} \\
    ${genome_fasta} \\
    ${reads_fwd} \\
    ${reads_rev} | \\
    \\
    sambamba view \\
      --sam-input \\
      --format bam \\
      --compression-level 0 \\
      --nthreads ${task.cpus} \\
      /dev/stdin | \\
    \\
    sambamba sort \\
      --nthreads ${task.cpus} \\
      --out ${meta.split}.${meta.sample_id}.${meta.read_group}.bam \\
      /dev/stdin
    """

    // TODO(SW): Versions.
    // """
    // touch bar

    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     bwamem2: foo
    // END_VERSIONS
    // """

    stub:
    """
    touch ${meta.split}.${meta.sample_id}.${meta.read_group}.bam
    """

    // TODO(MV): Versions.
    // """
    // touch bar
    // echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    // """
}
