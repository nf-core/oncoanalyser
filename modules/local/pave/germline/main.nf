process PAVE_GERMLINE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-pave:1.8--hdfd78af_1' :
        'biocontainers/hmftools-pave:1.8--hdfd78af_1' }"

    input:
    tuple val(meta), path(sage_vcf), path(sage_tbi)
    path genome_fasta
    val genome_ver
    path genome_fai
    path sage_blocklist_regions
    path sage_blocklist_sites
    path clinvar_annotations
    path segment_mappability
    path driver_gene_panel
    path ensembl_data_resources

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: index
    path 'versions.yml'                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    pave \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -input_vcf ${sage_vcf} \\
        -output_vcf ${meta.sample_id}.pave.germline.vcf.gz \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -clinvar_vcf ${clinvar_annotations} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -mappability_bed ${segment_mappability} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -blacklist_bed ${sage_blocklist_regions} \\
        -blacklist_vcf ${sage_blocklist_sites} \\
        -gnomad_no_filter \\
        -threads ${task.cpus} \\
        -output_dir ./ \\
        -log_level ${params.module_log_level}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pave: \$(pave -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.pave.germline.vcf.gz{,.tbi}

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
