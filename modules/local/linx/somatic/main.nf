process LINX_SOMATIC {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/linx:1.24.1--0'

    input:
    tuple val(meta), path(purple_dir)
    val genome_ver
    path ensembl_data_resources
    path known_fusion_data
    path driver_gene_panel
    path gene_id_file

    output:
    tuple val(meta), path('linx_somatic/'), emit: annotation_dir
    path 'versions.yml'                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def gene_id_file_arg = gene_id_file ? "-gene_id_file ${gene_id_file}" : ''

    """
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -sample ${meta.id} \\
            -sv_vcf ${purple_dir}/${meta.id}.purple.sv.vcf.gz \\
            -purple_dir ${purple_dir} \\
            ${gene_id_file_arg} \\
            -ref_genome_version ${genome_ver} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -known_fusion_file ${known_fusion_data} \\
            -driver_gene_panel ${driver_gene_panel} \\
            -write_vis_data \\
            -output_dir linx_somatic/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        linx: \$(java -jar ${task.ext.jarPath} | sed 's/^.*Linx version: //')
    END_VERSIONS
    """

    stub:
    """
    mkdir linx_somatic/
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
