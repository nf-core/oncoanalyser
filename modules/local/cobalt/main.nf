process COBALT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-cobalt:2.0--hdfd78af_0' :
        'biocontainers/hmftools-cobalt:2.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai)
    path gc_profile
    path diploid_regions
    path target_region_normalisation

    output:
    tuple val(meta), path('cobalt/'), emit: cobalt_dir
    path 'versions.yml'             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def reference_arg = meta.containsKey('normal_id') ? "-reference ${meta.normal_id}" : ''
    def reference_bam_arg = normal_bam ? "-reference_bam ${normal_bam}" : ''

    def diploid_regions_arg = diploid_regions ? "-tumor_only_diploid_bed ${diploid_regions}" : ''
    def target_region_arg = target_region_normalisation ? "-target_region ${target_region_normalisation}" : ''

    def run_mode = Utils.getEnumFromString(params.mode, Constants.RunMode)
    def pcf_gamma_arg = run_mode === Constants.RunMode.TARGETED && !meta.containsKey('normal_id')
        ? "-pcf_gamma 50" : ""

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
        ${target_region_arg} \\
        ${pcf_gamma_arg} \\
        -output_dir cobalt/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobalt: \$(cobalt -version | sed -n '/^Cobalt version/ { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p cobalt/
    touch cobalt/placeholder

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
