process SAMTOOLS_MERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(alns)
    path genome_fasta
    val format

    output:
    tuple val(meta), path("${meta.sample_id}.plain.*")        , topic: plain_aln
    tuple val(meta), val('samtools_merge'), path('.command.*'), topic: command_files
    path 'versions.yml'                                       , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def format_args = ''
    if (format == 'cram') {
        format_args = "-O cram --output-fmt-option version=3.0 --output-fmt-option reference=${genome_fasta} --output-fmt-option store_nm=1"
    }

    def format_index
    if (format == 'cram') {
        format_index = 'crai'
    } else if (format == 'bam') {
        format_index = 'bai'
    } else {
        error "did not receive either 'cram' or 'bam' for the 'format' input value"
    }

    """
    samtools merge \\
        ${format_args} \\
        --write-index \\
        --threads ${task.cpus} \\
        -o "${meta.sample_id}.plain.${format}##idx##${meta.sample_id}.plain.${format}.${format_index}" \\
        ${alns.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def output_fn = "${meta.sample_id}.plain"
    def index_ext = format == 'cram' ? 'crai' : 'bai'

    """
    touch ${output_fn}.${format}
    touch ${output_fn}.${format}.${index_ext}

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
