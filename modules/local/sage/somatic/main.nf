// NOTE(SW): logic that determines BQR outputs assumes '-output_vcf' is a path that includes at least leading one directory

process SAGE_SOMATIC {
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
    path sage_known_hotspots_somatic
    path sage_actionable_panel
    path sage_coverage_panel
    path sage_highconf_regions
    path ensembl_data_resources

    output:
    tuple val(meta), path('somatic/*.sage.somatic.vcf.gz'), path('somatic/*.sage.somatic.vcf.gz.tbi'), emit: vcf
    tuple val(meta), path('somatic/')                                                                , emit: sage_dir
    path 'versions.yml'                                                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def reference_arg = meta.containsKey('normal_id') ? "-reference ${meta.normal_id}" : ''
    def reference_bam_arg = normal_bam ? "-reference_bam ${normal_bam}" : ''

    """
    mkdir -p somatic/

    sage \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
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
        -output_vcf somatic/${meta.tumor_id}.sage.somatic.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: \$(sage -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p somatic/
    touch somatic/${meta.tumor_id}.sage.somatic.vcf.gz
    touch somatic/${meta.tumor_id}.sage.somatic.vcf.gz.tbi
    touch somatic/${meta.tumor_id}.gene.coverage.tsv
    touch somatic/${meta.tumor_id}.sage.bqr.png
    touch somatic/${meta.tumor_id}.sage.bqr.tsv
    touch somatic/${meta.normal_id}.sage.bqr.png
    touch somatic/${meta.normal_id}.sage.bqr.tsv
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
