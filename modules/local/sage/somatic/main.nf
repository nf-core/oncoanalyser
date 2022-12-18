// NOTE(SW): logic that determines BQR outputs assumes '-out' is a path that includes at least one directory

process SAGE_SOMATIC {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/sage:3.0.3--0'

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai)
    path genome_fasta
    path genome_fai
    path genome_dict
    val genome_ver
    path sage_known_hotspots_somatic
    path sage_coding_panel
    path sage_highconf_regions
    path ensembl_data_resources

    output:
    tuple val(meta), path("${meta.tumor_id}.sage_somatic.vcf.gz"), emit: vcf
    path '*gene.coverage.tsv'                                    , emit: gene_coverage, optional: true
    path '*sage.bqr.png'                                         , emit: bqr_png, optional: true
    path '*sage.bqr.tsv'                                         , emit: bqr_tsv, optional: true
    path 'versions.yml'                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    java \\
        -Xmx${task.memory.giga}g \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -reference ${meta.normal_id} \\
            -reference_bam ${normal_bam} \\
            -tumor ${meta.tumor_id} \\
            -tumor_bam ${tumor_bam} \\
            -ref_genome_version ${genome_ver} \\
            -ref_genome ${genome_fasta} \\
            -hotspots ${sage_known_hotspots_somatic} \\
            -panel_bed ${sage_coding_panel} \\
            -coverage_bed ${sage_coding_panel} \\
            -high_confidence_bed ${sage_highconf_regions} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -write_bqr_data \\
            -write_bqr_plot \\
            -threads ${task.cpus} \\
            -out ./${meta.tumor_id}.sage_somatic.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: \$(java -jar ${task.ext.jarPath} | head -n1 | sed 's/.*Sage version: //')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.tumor_id}.sage_somatic.vcf.gz"
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
