process SAGE_GERMLINE {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/sage:3.0.3--0'

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai)
    path genome_fasta
    path genome_fai
    path genome_dict
    val genome_ver
    path sage_known_hotspots_germline
    path sage_coding_panel
    path sage_highconf_regions
    path ensembl_data_resources

    output:
    tuple val(meta), path("${meta.tumor_id}.sage_germline.vcf.gz"), emit: vcf
    path 'versions.yml'                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    java \\
        -Xmx${task.memory.giga}g \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -reference ${meta.tumor_id} \\
            -reference_bam ${tumor_bam} \\
            -tumor ${meta.normal_id} \\
            -tumor_bam ${normal_bam} \\
            -ref_genome_version ${genome_ver} \\
            -ref_genome ${genome_fasta} \\
            -hotspots ${sage_known_hotspots_germline} \\
            -panel_bed ${sage_coding_panel} \\
            -high_confidence_bed ${sage_highconf_regions} \\
            -ensembl_data_resources ${ensembl_data_resources} \\
            -hotspot_min_tumor_qual 50 \\
            -panel_min_tumor_qual 75 \\
            -hotspot_max_germline_vaf 100 \\
            -hotspot_max_germline_rel_raw_base_qual 100 \\
            -panel_max_germline_vaf 100 \\
            -panel_max_germline_rel_raw_base_qual 100 \\
            -mnv_filter_enabled false \\
            -panel_only \\
            -threads ${task.cpus} \\
            -out ${meta.tumor_id}.sage_germline.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: \$(java -jar ${task.ext.jarPath} | head -n1 | sed 's/.*Sage version: //')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.tumor_id}.sage_germline.vcf.gz"
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
