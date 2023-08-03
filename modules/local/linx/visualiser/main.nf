process LINX_VISUALISER {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/linx:1.24.1--0'

    input:
    tuple val(meta), path(linx_annotation_dir)
    val genome_ver
    path ensembl_data_resources

    output:
    tuple val(meta), path('linx_visualiser/plot/'), emit: visualiser_dir
    path 'versions.yml'                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    # NOTE(SW): LINX v1.23.2 and current latest version (v1.23.6) require trailing slashes for the -plot_out and -data_out arguments since no filesystem separator is used when constructing fusion plot output filepaths.
    # https://github.com/hartwigmedical/hmftools/blob/linx-v1.23.6/linx/src/main/java/com/hartwig/hmftools/linx/visualiser/circos/ChromosomeRangeExecution.java#L22-L29
    # https://github.com/hartwigmedical/hmftools/blob/linx-v1.23.6/linx/src/main/java/com/hartwig/hmftools/linx/visualiser/circos/FusionExecution.java#L18-L23

    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -cp ${task.ext.jarPath} \\
        com.hartwig.hmftools.linx.visualiser.SvVisualiser \\
            ${args} \\
            -sample ${meta.id} \\
            -vis_file_dir ${linx_annotation_dir} \\
            -ref_genome_version ${genome_ver} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -circos ${task.ext.circosPath} \\
            -threads ${task.cpus} \\
            -plot_out linx_visualiser/plot/ \\
            -data_out linx_visualiser/data/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        linx: \$(java -jar ${task.ext.jarPath} | sed 's/^.*Linx version: //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p linx_visualiser/plot/
    echo -e '${task.process}:\n  stub: noversions\n' > versions.yml
    """
}
