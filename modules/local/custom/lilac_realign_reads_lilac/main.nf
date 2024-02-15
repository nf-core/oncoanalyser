process CUSTOM_REALIGNREADS {
    tag "${meta.id}"
    label 'process_low'

    conda "bwa=0.7.17 samtools=1.19.2 sambamba=1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4908f45b5b676e6a2fed6b1977d445b16b7b8dee:2411d64e0a784487e81828123aaf68a549531e5c-0' :
        'quay.io/biocontainers/mulled-v2-4908f45b5b676e6a2fed6b1977d445b16b7b8dee:2411d64e0a784487e81828123aaf68a549531e5c-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path reference
    path reference_indices

    output:
    tuple val(meta), path("*realigned.bam"), path("*realigned.bam.bai"), emit: bam
    path 'versions.yml'                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    sambamba sort -n ${bam} -o ${meta.sample_id}_sorted.bam
    samtools fastq -@${task.threads} ${meta.sample_id}_sorted.bam \\
            -1 ${meta.sample_id}_R1.fastq.gz \\
            -2 ${meta.sample_id}_R2.fastq.gz \\
            -0 ${meta.sample_id}_other.fastq.gz \\
            -s ${meta.sample_id}_singleton.fastq.gz;

    bwa mem \\
        -t${task.cpus} \\
        -Y \\
        ${reference} \\
        ${meta.sample_id}_R1.fastq.gz \\
        ${meta.sample_id}_R2.fastq.gz | \\
        samtools sort -T tmp -o ${bam.baseName}.realigned.bam
    samtools index ${bam.baseName}.realigned.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        sambamba: \$(sambamba --version 2>&1 | sed -n '/sambamba/ s/^sambamba \\(.\\+\\)/\\1/p' | head -n1)
    END_VERSIONS
    """

    stub:
    """
    touch ${bam.baseName}.realigned.bam ${bam.baseName}.realigned.bam.bai
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
