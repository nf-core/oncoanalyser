process SAMBAMBA_INDEX {
    tag "${meta.subject_id}__${meta.sample_id}"

    container 'docker.io/scwatts/sambamba:1.0'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path('*bai'), emit: bam
    path 'versions.yml'                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    sambamba index \\
      --nthreads ${task.cpus} \\
      ${meta.split}.${meta.sample_id}.${meta.read_group}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sambamba: \$(sambamba --version 2>&1 | egrep '^sambamba' | head -n 1 | awk '{ print \$NF }')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.split}.${meta.sample_id}.${meta.read_group}.bam.bai

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
