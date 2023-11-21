process LINX_VISUALISER {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/linx:1.24.1--0'

    input:
    tuple val(meta), path(linx_annotation_dir)
    val genome_ver
    path ensembl_data_resources

    output:
    tuple val(meta), path('linx_visualiser/all/')       , emit: visualiser_dir_all
    tuple val(meta), path('linx_visualiser/reportable/'), emit: visualiser_dir_reportable
    path 'versions.yml'                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    # NOTE(SW): LINX v1.24.1 require trailing slashes for the -plot_out and -data_out arguments since no filesystem
    # separator is used when constructing fusion plot output filepaths.

    # https://github.com/hartwigmedical/hmftools/blob/linx-v1.24.1/linx/src/main/java/com/hartwig/hmftools/linx/visualiser/circos/ChromosomeRangeExecution.java#L22-L29
    # https://github.com/hartwigmedical/hmftools/blob/linx-v1.24.1/linx/src/main/java/com/hartwig/hmftools/linx/visualiser/circos/FusionExecution.java#L18-L23

    # Generate all chromosome and cluster plots by default

    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -cp ${task.ext.jarPath} \\
        com.hartwig.hmftools.linx.visualiser.SvVisualiser \\
            ${args} \\
            -sample ${meta.sample_id} \\
            -vis_file_dir ${linx_annotation_dir} \\
            -ref_genome_version ${genome_ver} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -circos ${task.ext.circosPath} \\
            -threads ${task.cpus} \\
            -plot_out linx_visualiser/all/ \\
            -data_out linx_visualiser/all_data/

    # Rerun LINX to render only reportable cluster plots in a separate directory. While this is regenerating existing
    # cluster plots, the number of reportable plots is generally very small and I prefer to rely on the internal LINX
    # logic to determine whether a cluster is reportable rather than attempting to infer manually to copy out target
    # plot files.

    # The ORANGE report receives only reportable clusters while the gpgr LINX report receives chromosome and all cluster
    # plots.

    # https://github.com/hartwigmedical/hmftools/blob/linx-v1.24.1/linx/src/main/java/com/hartwig/hmftools/linx/visualiser/SampleData.java#L220-L236

    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -cp ${task.ext.jarPath} \\
        com.hartwig.hmftools.linx.visualiser.SvVisualiser \\
            ${args} \\
            -sample ${meta.sample_id} \\
            -vis_file_dir ${linx_annotation_dir} \\
            -ref_genome_version ${genome_ver} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -circos ${task.ext.circosPath} \\
            -plot_reportable \\
            -threads ${task.cpus} \\
            -plot_out linx_visualiser/reportable/ \\
            -data_out linx_visualiser/reportable_data/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        linx: \$(java -jar ${task.ext.jarPath} | sed 's/^.*Linx version: //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p linx_visualiser/plot_{all,reportable}/
    touch linx_visualiser/plot_{all,reportable}/placeholder
    echo -e '${task.process}:\n  stub: noversions\n' > versions.yml
    """
}
