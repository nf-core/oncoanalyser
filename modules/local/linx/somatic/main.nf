process LINX_SOMATIC {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-linx:2.0.2--hdfd78af_0' :
        'biocontainers/hmftools-linx:2.0.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(purple_dir)
    val genome_ver
    path ensembl_data_resources
    path known_fusion_data
    path driver_gene_panel

    output:
    tuple val(meta), path('linx_somatic/'), emit: annotation_dir
    path 'versions.yml'                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    """
    linx \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -sv_vcf ${purple_dir}/${meta.sample_id}.purple.sv.vcf.gz \\
        -purple_dir ${purple_dir} \\
        -ref_genome_version ${genome_ver} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -known_fusion_file ${known_fusion_data} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -write_vis_data \\
        -write_neo_epitopes \\
        -output_dir linx_somatic/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        linx: \$(linx -version | sed -n '/^Linx version / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir linx_somatic/
    touch linx_somatic/placeholder

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
