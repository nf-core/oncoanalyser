process LINX_GERMLINE {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/linx:1.24.1--0'

    input:
    tuple val(meta), path(sv_vcf)
    val genome_ver
    path ensembl_data_resources
    path driver_gene_panel

    output:
    tuple val(meta), path('linx_germline/'), emit: annotation_dir
    path 'versions.yml'                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -sample ${meta.id} \\
            -sv_vcf ${sv_vcf} \\
            -germline \\
            -ref_genome_version ${genome_ver} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -driver_gene_panel ${driver_gene_panel} \\
            -output_dir linx_germline/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        linx: \$(java -jar ${task.ext.jarPath} | sed 's/^.*Linx version: //')
    END_VERSIONS
    """

    stub:
    """
    mkdir linx_germline/
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
