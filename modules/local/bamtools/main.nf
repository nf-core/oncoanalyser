process BAMTOOLS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-bam-tools:1.2.1--hdfd78af_0' :
        'quay.io/biocontainers/hmftools-bam-tools:1.2.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path genome_fasta
    val genome_ver

    output:
    tuple val(meta), path('*wgsmetrics'), emit: metrics
    path 'versions.yml'                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bamtools \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        com.hartwig.hmftools.bamtools.metrics.BamMetrics \\
        -sample ${meta.sample_id} \\
        -bam_file ${bam} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -threads ${task.cpus} \\
        -write_old_style \\
        -log_level INFO \\
        -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamtools: \$(bamtools -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.wgsmetrics
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
