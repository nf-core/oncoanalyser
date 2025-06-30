process TEAL_PIPELINE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-teal:1.3.5--hdfd78af_0' :
        'biocontainers/hmftools-teal:1.3.5--hdfd78af_0' }"

    input:
    tuple val(meta),
        path(tumor_teal_bam), path(tumor_teal_bai),
        path(normal_teal_bam), path(normal_teal_bai),
        path(tumor_metrics_dir), path(normal_metrics_dir), path(cobalt_dir), path(purple_dir)
    val genome_ver

    output:
    tuple val(meta), path('teal/*.tsv*'), emit: teal_tsvs
    path 'versions.yml'                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def tumor_arg = tumor_teal_bam ? "-tumor ${meta.tumor_id}": ''
    def tumor_bam_arg = tumor_teal_bam ? "-tumor_bam ${tumor_teal_bam}": ''
    def tumor_wgs_metrics_arg = tumor_metrics_dir ? "-tumor_wgs_metrics ${tumor_metrics_dir}/${meta.tumor_id}.bam_metric.summary.tsv": ''
    def purple_arg = purple_dir ? "-purple ${purple_dir}": ''

    def reference_arg = normal_teal_bam ? "-reference ${meta.normal_id}" : ''
    def reference_bam_arg = normal_teal_bam ? "-reference_bam ${normal_teal_bam}" : ''
    def reference_wgs_metrics_arg = normal_metrics_dir ? "-reference_wgs_metrics ${normal_metrics_dir}/${meta.normal_id}.bam_metric.summary.tsv" : ''

    if (tumor_arg && ! purple_arg) error "TEAL requires PURPLE inputs when analysing tumor data"
    if (! tumor_arg && ! reference_arg) error "TEAL at least tumor or normal data for analyses"

    """
    teal \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        com.hartwig.hmftools.teal.TealPipelineApp \\
        ${args} \\
        ${reference_arg} \\
        ${reference_bam_arg} \\
        ${tumor_arg} \\
        ${tumor_bam_arg} \\
        -cobalt ${cobalt_dir} \\
        ${purple_arg} \\
        ${reference_wgs_metrics_arg} \\
        ${tumor_wgs_metrics_arg} \\
        -ref_genome_version ${genome_ver} \\
        -threads ${task.cpus} \\
        -output_dir teal/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        teal: \$(teal -version | sed -n '/Teal version/ { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p teal/

    ${ (meta.tumor_id != null) ? "touch teal/${meta.tumor_id}.teal.{tellength.tsv,breakend.tsv.gz}" : '' }
    ${ (meta.normal_id != null) ? "touch teal/${meta.normal_id}.teal.{tellength.tsv}" : '' }

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
