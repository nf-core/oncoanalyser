process QSEE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-qsee:1.0--hdfd78af_0' :
        'biocontainers/hmftools-qsee:1.0--hdfd78af_0' }"

    input:
    tuple val(meta),
        path(redux_tsvs_tumor, stageAs: "redux_tsvs_tumor/*"),
        path(redux_tsvs_normal, stageAs: "redux_tsvs_normal/*"),
        path(bamtools_dir_tumor, stageAs: 'bamtools_tumor'),
        path(bamtools_dir_normal, stageAs: 'bamtools_normal'),
        path(cobalt_dir),
        path(esvee_dir),
        path(purple_dir)
    path driver_gene_panel
    path cohort_percentiles
    val sequencing_platform
    val targeted_mode

    output:
    tuple val(meta), path('qsee/')                  , topic: qsee_dir
    tuple val(meta), val('qsee'), path('.command.*'), topic: command_files
    path 'versions.yml'                             , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def reference_arg = meta.containsKey('normal_id') ? "-reference ${meta.normal_id}" : ''
    def redux_ref_dir_arg = redux_tsvs_normal ? '-redux_ref_dir redux_tsvs_normal/' : ''
    def bamtools_ref_dir_arg = bamtools_dir_normal ? "-bam_metrics_ref_dir ${bamtools_dir_normal}" : ''

    def cobalt_dir_arg = cobalt_dir ? "-cobalt_dir ${cobalt_dir}" : ''
    def esvee_dir_arg = esvee_dir ? "-esvee_dir ${esvee_dir}" : ''

    def targeted_mode_arg = targeted_mode ? '-targeted_mode' : ''

    def cohort_percentiles_arg = ''
    if(! targeted_mode && sequencing_platform.toLowerCase() == 'illumina') {
        cohort_percentiles_arg = "-cohort_percentiles_file ${cohort_percentiles}"
    }

    """
    mkdir -p qsee/

    qsee \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -tumor ${meta.tumor_id} \\
        ${reference_arg} \\
        -redux_tumor_dir redux_tsvs_tumor/ \\
        ${redux_ref_dir_arg} \\
        -bam_metrics_tumor_dir ${bamtools_dir_tumor} \\
        ${bamtools_ref_dir_arg} \\
        ${cobalt_dir_arg} \\
        ${esvee_dir_arg} \\
        -purple_dir ${purple_dir} \\
        -sequencing_type ${sequencing_platform.toUpperCase()} \\
        -driver_gene_panel ${driver_gene_panel} \\
        ${cohort_percentiles_arg} \\
        ${targeted_mode_arg} \\
        -allow_missing_input \\
        ${log_level_arg} \\
        -output_dir \$(realpath qsee/)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qsee: \$(qsee -version | sed -n '/^Qsee version/ { s/^.* //p }')
        r: \$(R --version | sed -n '/^R version/ { s/^.*version //; s/ .*//p }')
        r-dplyr: \$(Rscript -e 'packageVersion("dplyr") |> as.character() |> writeLines()')
        r-ggplot2: \$(Rscript -e 'packageVersion("ggplot2") |> as.character() |> writeLines()')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p qsee/

    gzip <<< '' > qsee/${meta.tumor_id}.qsee.status.tsv.gz
    touch qsee/${meta.tumor_id}.qsee.vis.report.pdf
    gzip <<< '' > qsee/${meta.tumor_id}.qsee.vis.data.tsv.gz

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
