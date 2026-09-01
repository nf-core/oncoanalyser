process LILAC {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-lilac:2.0--hdfd78af_0' :
        'biocontainers/hmftools-lilac:2.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(normal_dna_aln), path(normal_dna_idx), path(tumor_dna_aln), path(tumor_dna_idx), path(tumor_rna_aln), path(tumor_rna_idx), path(purple_dir)
    path genome_fasta
    val genome_ver
    path genome_fai
    path lilac_resources, stageAs: 'lilac_resources'
    val targeted_mode
    val sequencing_platform

    output:
    tuple val(meta), path('lilac/')                  , topic: lilac_dir
    tuple val(meta), val('lilac'), path('.command.*'), topic: command_files
    path 'versions.yml'                              , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def sample_name = getSampleName(meta, tumor_dna_aln, normal_dna_aln)

    def normal_bam_arg = normal_dna_aln ? "-reference_bam ${normal_dna_aln}" : ''
    def tumor_dna_bam_arg = tumor_dna_aln ? "-tumor_bam ${tumor_dna_aln}" : ''
    def tumor_rna_bam_arg = tumor_rna_aln ? "-rna_bam ${tumor_rna_aln}" : ''

    def purple_dir_arg = purple_dir ? "-purple_dir ${purple_dir}" : ''

    def freq_score_penalty
    def gene_groups
    def targeted_arg

    if (targeted_mode) {
        freq_score_penalty = '0.0018'
        gene_groups = 'MHC_CLASS_1'
        targeted_arg = '-targeted_panel'
    } else {
        freq_score_penalty = '0.0009'
        gene_groups = 'ALL'
        targeted_arg = ''
    }

    """
    lilac \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${sample_name} \\
        ${normal_bam_arg} \\
        ${tumor_dna_bam_arg} \\
        ${tumor_rna_bam_arg} \\
        ${purple_dir_arg} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -resource_dir ${lilac_resources} \\
        -freq_score_penalty ${freq_score_penalty} \\
        -genes ${gene_groups} \\
        ${targeted_arg} \\
        -sequencing_type ${sequencing_platform.toUpperCase()} \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_dir lilac/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lilac: \$(lilac -version | sed -n '/^Lilac version / { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p lilac/

    touch lilac/.stub

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}

def getSampleName(meta, tumor_aln, normal_aln) {
    if (tumor_aln) {
        return meta.tumor_id
    } else if (normal_aln) {
        return meta.normal_id
    } else {
        error 'did not receive either a tumor or normal alignment'
    }
}
