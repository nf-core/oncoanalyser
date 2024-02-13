process BWA_MEM {
    tag "${meta.id}"

    // TODO(MC): What process label?
    // label 'process_medium'

    container 'docker.io/scwatts/bwa:0.7.17-sambamba'

    input:
    tuple val(meta), path(reads_fwd), path(reads_rev)
    path genome_fasta
    path genome_bwa_index

    output:
    tuple val(meta), path('*bam'), emit: bam

    // TODO(MC): How does this work?
    when:
    task.ext.when == null || task.ext.when

    // # TODO(MC): read group
    // # -R ${meta.read_group}

    script:
    """
    ln -s \$(find -L ${genome_bwa_index} -type f) ./

    bwa mem \\
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

    stub:
    """
    touch ${meta.split}.${meta.sample_id}.${meta.read_group}.bam
    """
}
