process CUSTOM_SLICE {
    tag "${meta.id}"
    label 'process_single'

    conda "samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h96c455f_1' :
        'biocontainers/samtools:1.21--h96c455f_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed

    output:
    tuple val(meta), path("*sliced.bam"), path("*sliced.bam.bai"), emit: bam
    path 'versions.yml'                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    samtools view \\
        ${args} \\
        --regions-file ${bed} \\
        -@${task.cpus} \\
        -Obam \\
        ${bam} | \\
        samtools sort -T tmp -o ${bam.baseName}.sliced.bam

    samtools index ${bam.baseName}.sliced.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed -n '/^samtools / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    touch ${bam.baseName}.sliced.bam ${bam.baseName}.sliced.bam.bai

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
