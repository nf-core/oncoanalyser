// NOTE(SW): use of tumor and normal sample names here is consistent with Pipeline5
//  - https://github.com/hartwigmedical/pipeline5/blob/v5.32/cluster/src/main/java/com/hartwig/pipeline/calling/sage/SageCommandBuilder.java#L95-L96
//  - https://github.com/hartwigmedical/pipeline5/blob/v5.32/cluster/src/main/java/com/hartwig/pipeline/calling/sage/SageCommandBuilder.java#L112-L118

// NOTE(SW): logic that determines BQR outputs assumes '-out' is a path that includes at least leading one directory

// TODO(SW): check whether intentional sample name switch for germline also affects BQR outputs

process SAGE_GERMLINE {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/sage:3.2.5--0'

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai)
    path genome_fasta
    val genome_ver
    path genome_fai
    path genome_dict
    path sage_known_hotspots_germline
    path sage_actionable_panel
    path sage_coverage_panel
    path sage_highconf_regions
    path ensembl_data_resources

    output:
    tuple val(meta), path('*.sage.germline.vcf.gz'), path('*.sage.germline.vcf.gz.tbi')                  , emit: vcf
    tuple val(meta), path('*.sage.germline.filtered.vcf.gz'), path('*.sage.germline.filtered.vcf.gz.tbi'), emit: vcf_filtered
    tuple val(meta), path("${meta.tumor_id}.sage.bqr.png")                                               , emit: tumor_bqr_png, optional: true
    tuple val(meta), path("${meta.normal_id}.sage.bqr.png")                                              , emit: normal_bqr_png, optional: true
    tuple val(meta), path('*gene.coverage.tsv')                                                          , emit: gene_coverage, optional: true
    path '*sage.bqr.tsv'                                                                                 , emit: bqr_tsv, optional: true
    path 'versions.yml'                                                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -tumor ${meta.normal_id} \\
            -tumor_bam ${normal_bam} \\
            -reference ${meta.tumor_id} \\
            -reference_bam ${tumor_bam} \\
            -ref_genome ${genome_fasta} \\
            -ref_genome_version ${genome_ver} \\
            -hotspots ${sage_known_hotspots_germline} \\
            -panel_bed ${sage_actionable_panel} \\
            -coverage_bed ${sage_coverage_panel} \\
            -high_confidence_bed ${sage_highconf_regions} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -hotspot_min_tumor_qual 50 \\
            -panel_min_tumor_qual 75 \\
            -hotspot_max_germline_vaf 100 \\
            -hotspot_max_germline_rel_raw_base_qual 100 \\
            -panel_max_germline_vaf 100 \\
            -panel_max_germline_rel_raw_base_qual 100 \\
            -ref_sample_count 0 \\
            -panel_only \\
            -write_bqr_data \\
            -write_bqr_plot \\
            -threads ${task.cpus} \\
            -out ./${meta.tumor_id}.sage.germline.vcf.gz

    bcftools view -f 'PASS' -o ${meta.tumor_id}.sage.germline.filtered.vcf.gz ${meta.tumor_id}.sage.germline.vcf.gz
    bcftools index -t ${meta.tumor_id}.sage.germline.filtered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: \$(java -jar ${task.ext.jarPath} | head -n1 | sed 's/.*Sage version: //')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.tumor_id}.sage.germline.vcf.gz"
    touch "${meta.tumor_id}.sage.germline.vcf.gz.tbi"
    touch "${meta.tumor_id}.sage.germline.filtered.vcf.gz"
    touch "${meta.tumor_id}.sage.germline.filtered.vcf.gz.tbi"
    touch "${meta.tumor_id}.sage.bqr.png"
    touch "${meta.tumor_id}.sage.bqr.tsv"
    touch "${meta.normal_id}.sage.bqr.png"
    touch "${meta.normal_id}.sage.bqr.tsv"
    touch "${meta.normal_id}.gene.coverage.tsv"
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
