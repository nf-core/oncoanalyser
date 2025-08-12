process COBALT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-cobalt:2.1--hdfd78af_1' :
        'biocontainers/hmftools-cobalt:2.1--hdfd78af_1' }"

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai)
    path gc_profile
    path diploid_regions
    path target_region_normalisation
    val is_targeted_mode

    output:
    tuple val(meta), path('cobalt/'), emit: cobalt_dir
    path 'versions.yml'             , emit: versions
    path '.command.*'               , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def reference_arg = meta.containsKey('normal_id') ? "-reference ${meta.normal_id}" : ''
    def reference_bam_arg = normal_bam ? "-reference_bam ${normal_bam}" : ''

    def target_region_norm_file_arg = target_region_normalisation ? "-target_region_norm_file ${target_region_normalisation}" : ''

    def is_tumor_only_mode = !meta.containsKey('normal_id')

    def pcf_gamma_arg = is_targeted_mode && is_tumor_only_mode
        ? "-pcf_gamma 50" : ""

    def diploid_regions_arg = !is_targeted_mode && is_tumor_only_mode
        ? "-tumor_only_diploid_bed ${diploid_regions}" : ""

    """
    cobalt \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -tumor ${meta.tumor_id} \\
        -tumor_bam ${tumor_bam} \\
        ${reference_arg} \\
        ${reference_bam_arg} \\
        -threads ${task.cpus} \\
        -gc_profile ${gc_profile} \\
        ${diploid_regions_arg} \\
        ${target_region_norm_file_arg} \\
        ${pcf_gamma_arg} \\
        ${log_level_arg} \\
        -output_dir cobalt/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobalt_run: \$(cobalt -version | sed -n '/^Cobalt version/ { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p cobalt/
    touch cobalt/placeholder

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
