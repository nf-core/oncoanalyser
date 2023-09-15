// NOTE(SW): logic that determines BQR outputs assumes '-out' is a path that includes at least leading one directory

process SAGE_SOMATIC {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/sage:3.3.1--0'

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai)
    path genome_fasta
    val genome_ver
    path genome_fai
    path genome_dict
    path sage_known_hotspots_somatic
    path sage_actionable_panel
    path sage_coverage_panel
    path sage_highconf_regions
    path ensembl_data_resources

    output:
    tuple val(meta), path('somatic/*.sage.somatic.vcf.gz'), path('somatic/*.sage.somatic.vcf.gz.tbi')                  , emit: vcf
    tuple val(meta), path('somatic/*.sage.somatic.filtered.vcf.gz'), path('somatic/*.sage.somatic.filtered.vcf.gz.tbi'), emit: vcf_filtered
    tuple val(meta), path('somatic/')                                                                                  , emit: sage_dir
    path 'versions.yml'                                                                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def reference_arg = meta.containsKey('normal_id') ? "-reference ${meta.normal_id}" : ''
    def reference_bam_arg = normal_bam ? "-reference_bam ${normal_bam}" : ''

    """
    mkdir -p somatic/

    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            ${reference_arg} \\
            ${reference_bam_arg} \\
            -tumor ${meta.tumor_id} \\
            -tumor_bam ${tumor_bam} \\
            -ref_genome ${genome_fasta} \\
            -ref_genome_version ${genome_ver} \\
            -hotspots ${sage_known_hotspots_somatic} \\
            -panel_bed ${sage_actionable_panel} \\
            -coverage_bed ${sage_coverage_panel} \\
            -high_confidence_bed ${sage_highconf_regions} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -write_bqr_data \\
            -write_bqr_plot \\
            -threads ${task.cpus} \\
            -out somatic/${meta.tumor_id}.sage.somatic.vcf.gz

    bcftools view \\
        -f 'PASS' \\
        -o somatic/${meta.tumor_id}.sage.somatic.filtered.vcf.gz \\
        somatic/${meta.tumor_id}.sage.somatic.vcf.gz

    bcftools index -t somatic/${meta.tumor_id}.sage.somatic.filtered.vcf.gz

    # NOTE(SW): hard coded since there is no reliable way to obtain version information.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: 3.3
    END_VERSIONS
    """

    stub:
    """
    mkdir -p somatic/
    touch somatic/${meta.tumor_id}.sage.somatic.vcf.gz
    touch somatic/${meta.tumor_id}.sage.somatic.vcf.gz.tbi
    touch somatic/${meta.tumor_id}.sage.somatic.filtered.vcf.gz
    touch somatic/${meta.tumor_id}.sage.somatic.filtered.vcf.gz.tbi
    touch somatic/${meta.tumor_id}.gene.coverage.tsv
    touch somatic/${meta.tumor_id}.sage.bqr.png
    touch somatic/${meta.tumor_id}.sage.bqr.tsv
    touch somatic/${meta.normal_id}.sage.bqr.png
    touch somatic/${meta.normal_id}.sage.bqr.tsv
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
