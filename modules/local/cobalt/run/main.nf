process COBALT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-cobalt:3.0--hdfd78af_0' :
        'biocontainers/hmftools-cobalt:3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_aln), path(tumor_idx), path(normal_aln), path(normal_idx)
    path genome_fasta
    val genome_ver
    path genome_fai
    path gc_profile
    path diploid_regions
    path target_regions_normalisation
    val targeted_mode

    output:
    tuple val(meta), path('cobalt/')                  , topic: cobalt_dir
    tuple val(meta), val('cobalt'), path('.command.*'), topic: command_files
    path 'versions.yml'                               , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def reference_arg = meta.containsKey('normal_id') ? "-reference ${meta.normal_id}" : ''
    def reference_bam_arg = normal_aln ? "-reference_bam ${normal_aln}" : ''

    def target_regions_norm_file_arg = target_regions_normalisation ? "-target_region_norm_file ${target_regions_normalisation}" : ''

    def tumor_only_mode = ! meta.containsKey('normal_id')

    def pcf_gamma_arg = targeted_mode && tumor_only_mode ? '-pcf_gamma 50' : ''

    def diploid_regions_arg = ! targeted_mode && tumor_only_mode ? "-tumor_only_diploid_bed ${diploid_regions}" : ''

    """
    cobalt \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -tumor ${meta.tumor_id} \\
        -tumor_bam ${tumor_aln} \\
        ${pcf_gamma_arg} \\
        ${reference_arg} \\
        ${reference_bam_arg} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -gc_profile ${gc_profile} \\
        ${diploid_regions_arg} \\
        ${target_regions_norm_file_arg} \\
        ${log_level_arg} \\
        -threads ${task.cpus} \\
        -output_dir cobalt/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobalt: \$(cobalt -version | sed -n '/^Cobalt version/ { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
        r: \$(R --version | sed -n '/^R version/ { s/^.*version //; s/ .*//p }')
        r-dplyr: \$(Rscript -e 'packageVersion("dplyr") |> as.character() |> writeLines()')
        bioconductor-copynumber: \$(Rscript -e 'packageVersion("copynumber") |> as.character() |> writeLines()')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p cobalt/

    touch cobalt/.stub

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
