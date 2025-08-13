process NEO_ANNOTATE_FUSIONS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-isofox:1.7.2--hdfd78af_1' :
        'biocontainers/hmftools-isofox:1.7.2--hdfd78af_1' }"

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
    path '.command.*'                           , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    """
    mkdir -p isofox/

    isofox \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -bam_file ${bam} \\
        -functions NEO_EPITOPES \\
        -read_length ${read_length} \\
        -neo_dir ${neo_finder_dir} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isofox: \$(isofox -version | sed -n '/^Isofox version / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.isf.neoepitope.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
