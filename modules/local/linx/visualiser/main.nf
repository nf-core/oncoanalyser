process LINX_VISUALISER {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/linx:1.23.2--0'

    input:
    tuple val(meta), path(linx)
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
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -cp ${task.ext.jarPath} \\
        com.hartwig.hmftools.linx.visualiser.SvVisualiser \\
            ${args} \\
            -sample ${meta.id} \\
            -vis_file_dir ${linx} \\
            -ref_genome_version ${genome_ver} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -circos ${task.ext.path_circos} \\
            -threads ${task.cpus} \\
            -plot_out linx_visualiser/plot \\
            -data_out linx_visualiser/data

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        linx: \$(java -jar ${task.ext.jarPath} | sed 's/^.*LINX version: //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p linx_visualiser/plot/
    echo -e '${task.process}:\n  stub: noversions\n' > versions.yml
    """
}
