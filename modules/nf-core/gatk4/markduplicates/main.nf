nextflow.enable.types = true

process GATK4_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.6.1.0--py310hdfd78af_0' :
        'biocontainers/gatk4:4.6.1.0--py310hdfd78af_0' }"

    input:
    tuple(meta: Map, bam: Path)
    fasta: Path?
    fasta_fai: Path?

    topic:
    tuple(meta, file('*.cram', optional: true))              >> 'gatk4_markduplicates_cram'
    tuple(meta, file('*.bam', optional: true))               >> 'gatk4_markduplicates_bam'
    tuple(meta, file('*.crai', optional: true))              >> 'gatk4_markduplicates_crai'
    tuple(meta, file('*.bai', optional: true))               >> 'gatk4_markduplicates_bai'
    tuple(meta, file('*.metrics'))                           >> 'gatk4_markduplicates_metrics'
    tuple(meta, 'gatk4_markduplicates', files('.command.*')) >> 'command_files'
    file('versions.yml')                                     >> 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    prefix = task.ext.prefix ?: "${meta.sample_id}"

    def input_list = bam.collect { s -> "--INPUT ${s}" }.join(' ')
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M" MarkDuplicates \\
        $input_list \\
        --OUTPUT ${prefix}.md.bam \\
        --METRICS_FILE ${prefix}.md.metrics \\
        --TMP_DIR . \\
        --CREATE_INDEX \\
        ${reference} \\
        $args

    mv ${prefix}.md.bai ${prefix}.md.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.sample_id}"

    """
    touch ${prefix}.md.bam
    touch ${prefix}.md.bam.bai
    touch ${prefix}.md.metrics

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
