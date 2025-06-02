process ISOFOX {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-isofox:1.7.1--hdfd78af_1' :
        'biocontainers/hmftools-isofox:1.7.1--hdfd78af_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val functions
    val read_length
    path genome_fasta
    val genome_ver
    path genome_fai
    path ensembl_data_resources
    path known_fusion_data
    path exp_counts
    path exp_gc_ratios
    path gene_ids
    path tpm_norm

    output:
    tuple val(meta), path('isofox/'), emit: isofox_dir
    path 'versions.yml'             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def functions_arg = functions ? "-functions \'${functions}\'" : ''

    def exp_counts_arg = exp_counts ? "-exp_counts_file ${exp_counts}" : ''
    def exp_gc_ratios_arg = exp_gc_ratios ? "-exp_gc_ratios_file ${exp_gc_ratios}" : ''

    def gene_ids_arg = gene_ids ? "-gene_id_file ${gene_ids}" : ''
    def tpm_norm_arg = tpm_norm ? "-panel_tpm_norm_file ${tpm_norm}" : ''

    """
    mkdir -p isofox/

    isofox \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -bam_file ${bam} \\
        ${functions_arg} \\
        -read_length ${read_length} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -known_fusion_file ${known_fusion_data} \\
        ${exp_counts_arg} \\
        ${exp_gc_ratios_arg} \\
        ${gene_ids_arg} \\
        ${tpm_norm_arg} \\
        -threads ${task.cpus} \\
        -output_dir isofox/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isofox: \$(isofox -version | sed -n '/^Isofox version / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p isofox/
    touch isofox/placeholder

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
