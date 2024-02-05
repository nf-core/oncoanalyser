process ANNOTATE_FUSIONS {
    tag "${meta.id}"
    label 'process_medium'

    container 'quay.io/biocontainers/hmftools-isofox:1.7.1--hdfd78af_0'

    input:
    tuple val(meta), path(neo_finder_dir), path(bam), path(bai)
    val read_length
    path genome_fasta
    val genome_ver
    path genome_fai
    path ensembl_data_resources

    output:
    tuple val(meta), path('*isf.neoepitope.tsv'), emit: annotated_fusions
    path 'versions.yml'                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p isofox/

    isofox \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -bam_file ${bam} \\
        -functions NEO_EPITOPES \\
        -neo_dir ${neo_finder_dir} \\
        -read_length ${read_length} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -threads ${task.cpus} \\
        -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isofox: \$(isofox -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.isf.neoepitope.tsv
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
