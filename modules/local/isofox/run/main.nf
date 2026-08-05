process ISOFOX {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-isofox:2.0.1--hdfd78af_0' :
        'biocontainers/hmftools-isofox:2.0.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(aln), path(idx)
    val functions
    val read_length
    path genome_fasta
    val genome_ver
    path genome_fai
    path ensembl_data_resources
    path driver_gene_panel
    path known_fusion_data
    path excluded_regions
    path gene_distribution
    path alt_sj_distribution
    path exp_counts
    path exp_gc_ratios
    path tpm_norm

    output:
    tuple val(meta), path('isofox/')                  , topic: isofox_dir
    tuple val(meta), val('isofox'), path('.command.*'), topic: command_files
    path 'versions.yml'                               , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def functions_arg = functions ? "-functions \'${functions}\'" : ''

    def exp_counts_arg = exp_counts ? "-exp_counts_file ${exp_counts}" : ''
    def exp_gc_ratios_arg = exp_gc_ratios ? "-exp_gc_ratios_file ${exp_gc_ratios}" : ''

    def tpm_norm_arg = tpm_norm ? "-panel_tpm_norm_file ${tpm_norm}" : ''

    """
    mkdir -p isofox/

    isofox \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        ${functions_arg} \\
        -read_length ${read_length} \\
        -bam_file ${aln} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -known_fusion_file ${known_fusion_data} \\
        -excluded_regions ${excluded_regions} \\
        -gene_distribution_file ${gene_distribution} \\
        -alt_sj_cohort_file ${alt_sj_distribution} \\
        ${exp_counts_arg} \\
        ${exp_gc_ratios_arg} \\
        ${tpm_norm_arg} \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_dir isofox/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isofox: \$(isofox -version | sed -n '/^Isofox version / { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p isofox/

    touch isofox/.stub

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
