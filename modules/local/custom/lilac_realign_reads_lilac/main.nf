process CUSTOM_REALIGNREADS {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.21 bioconda::sambamba=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4dde50190ae599f2bb2027cb2c8763ea00fb5084:b39930b2feee9b7a3e5cbef88b34f0d5c5f64897-0' :
        'biocontainers/mulled-v2-4dde50190ae599f2bb2027cb2c8763ea00fb5084:b39930b2feee9b7a3e5cbef88b34f0d5c5f64897-0' }"

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
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def args5 = task.ext.args5 ?: ''

    """
    sambamba sort ${args} -n ${bam} -o ${meta.sample_id}_sorted.bam

    samtools fastq ${args2} -@${task.threads} ${meta.sample_id}_sorted.bam \\
            -1 ${meta.sample_id}_R1.fastq.gz \\
            -2 ${meta.sample_id}_R2.fastq.gz \\
            -0 ${meta.sample_id}_other.fastq.gz \\
            -s ${meta.sample_id}_singleton.fastq.gz;

    bwa-mem2 mem \\
        ${args3} \\
        -Y \\
        -t ${task.cpus} \\
        ${reference} \\
        ${meta.sample_id}_R1.fastq.gz \\
        ${meta.sample_id}_R2.fastq.gz | \\
        \\
        sambamba view \\
            ${args4} \\
            --sam-input \\
            --format bam \\
            --compression-level 0 \\
            --nthreads ${task.cpus} \\
            /dev/stdin | \\
        \\
        sambamba sort \\
            ${args5} \\
            --nthreads ${task.cpus} \\
            --out ${bam.baseName}.realigned.bam \\
            /dev/stdin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version)
        sambamba: \$(sambamba --version 2>&1 | sed -n '/^sambamba / { s/^.* //p }' | head -n1)
        samtools: \$(samtools --version | sed -n '/^samtools / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    touch ${bam.baseName}.realigned.bam ${bam.baseName}.realigned.bam.bai

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
