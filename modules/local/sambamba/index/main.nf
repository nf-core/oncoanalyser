process SAMBAMBA_INDEX {
    tag "${meta.id}"

    container 'docker.io/scwatts/sambamba:1.0'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path('*bai'), emit: bam

    script:
    """
    sambamba index \\
      --nthreads ${task.cpus} \\
      ${meta.split}.${meta.sample_id}.${meta.read_group}.bam
    """

    stub:
    """
    touch ${meta.split}.${meta.sample_id}.${meta.read_group}.bam.bai
    """
}
