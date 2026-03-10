process QSEE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-qsee:1.0--hdfd78af_0' :
        'biocontainers/hmftools-qsee:1.0--hdfd78af_0' }"

    input:
    tuple val(meta),
        path(redux_somatic_tsv, stageAs: "redux_somatic/*"),
        path(redux_germline_tsv, stageAs: "redux_germline/*"),
        path(bamtools_somatic_dir, stageAs: 'bamtools_somatic'),
        path(bamtools_germline_dir, stageAs: 'bamtools_germline'),
        path(cobalt_dir),
        path(esvee_dir),
        path(purple_dir)
    path driver_gene_panel
    path cohort_percentiles
    val targeted_mode

    output:
    tuple val(meta), path('qsee/'), emit: qsee_dir
    path 'versions.yml'           , emit: versions
    path '.command.*'             , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def reference_arg = meta.normal_id ? "-reference ${meta.normal_id}" : ''
    def redux_ref_dir_arg = redux_germline_tsv ? "-redux_ref_dir redux_germline/" : ''
    def bamtools_ref_dir_arg = bamtools_germline_dir ? "-ref_metrics_dir ${bamtools_germline_dir}" : ''

    def cobalt_dir_arg = cobalt_dir ? "-cobalt_dir ${cobalt_dir}" : ''
    def esvee_dir_arg = esvee_dir ? "-esvee_dir ${esvee_dir}" : ''

    def cohort_percentiles_arg = !targeted_mode ? "-cohort_percentiles_file ${cohort_percentiles}" : ''

    """
    mkdir -p qsee/

    qsee \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -tumor ${meta.tumor_id} \\
        ${reference_arg} \\
        -redux_tumor_dir redux_somatic/ \\
        ${redux_ref_dir_arg} \\
        -tumor_metrics_dir ${bamtools_somatic_dir} \\
        ${bamtools_ref_dir_arg} \\
        ${cobalt_dir_arg} \\
        ${esvee_dir_arg} \\
        -purple_dir ${purple_dir} \\
        -driver_gene_panel ${driver_gene_panel} \\
        ${cohort_percentiles_arg} \\
        -allow_missing_input \\
        ${log_level_arg} \\
        -output_dir \$(realpath qsee/)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qsee: \$(qsee -version | sed -n '/^Qsee version/ { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p qsee/

    touch qsee/${meta.tumor_id}.qsee.status.tsv.gz
    touch qsee/${meta.tumor_id}.qsee.vis.report.pdf
    touch qsee/${meta.tumor_id}.qsee.vis.data.tsv.gz

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
