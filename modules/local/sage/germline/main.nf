process SAGE_GERMLINE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sage:4.1--hdfd78af_0' :
        'biocontainers/hmftools-sage:4.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai), path(redux_tsvs)
    path genome_fasta
    val genome_ver
    path genome_fai
    path genome_dict
    path sage_known_hotspots_germline
    path sage_highconf_regions
    path driver_gene_panel
    path ensembl_data_resources
    val targeted_mode

    output:
    tuple val(meta), path('germline/*.sage.germline.vcf.gz'), path('germline/*.sage.germline.vcf.gz.tbi'), emit: vcf
    tuple val(meta), path('germline/')                                                                   , emit: sage_dir
    path 'versions.yml'                                                                                  , emit: versions
    path '.command.*'                                                                                    , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def high_depth_mode_arg = targeted_mode ? '-high_depth_mode' : ''

    """
    mkdir -p germline/

    sage \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        -tumor ${meta.normal_id} \\
        -tumor_bam ${normal_bam} \\
        -reference ${meta.tumor_id} \\
        -reference_bam ${tumor_bam} \\
        -jitter_param_dir ./ \\
        -ref_sample_count 0 \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -hotspots ${sage_known_hotspots_germline} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -high_confidence_bed ${sage_highconf_regions} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -germline \\
        -panel_only \\
        ${high_depth_mode_arg} \\
        -bqr_write_plot \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
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
