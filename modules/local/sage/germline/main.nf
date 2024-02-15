process SAGE_GERMLINE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sage:3.4--hdfd78af_1' :
        'quay.io/biocontainers/hmftools-sage:3.4--hdfd78af_1' }"

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
    tuple val(meta), path('germline/*.sage.germline.vcf.gz'), path('germline/*.sage.germline.vcf.gz.tbi'), emit: vcf
    tuple val(meta), path('germline/')                                                                   , emit: sage_dir
    path 'versions.yml'                                                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir germline/

    sage \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
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
        -output_vcf germline/${meta.tumor_id}.sage.germline.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: \$(sage -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p germline/
    touch germline/${meta.tumor_id}.sage.germline.vcf.gz
    touch germline/${meta.tumor_id}.sage.germline.vcf.gz.tbi
    touch germline/${meta.tumor_id}.sage.bqr.png
    touch germline/${meta.tumor_id}.sage.bqr.tsv
    touch germline/${meta.normal_id}.sage.bqr.png
    touch germline/${meta.normal_id}.sage.bqr.tsv
    touch germline/${meta.normal_id}.gene.coverage.tsv
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
