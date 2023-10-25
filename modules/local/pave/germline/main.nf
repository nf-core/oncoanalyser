// NOTE(SW): use of tumor sample name here is consistent with Pipeline5
//  - https://github.com/hartwigmedical/pipeline5/blob/v5.32/cluster/src/main/java/com/hartwig/pipeline/tertiary/pave/PaveGermline.java#L35-L39
//  - https://github.com/hartwigmedical/pipeline5/blob/v5.32/cluster/src/main/java/com/hartwig/pipeline/tertiary/pave/PaveArguments.java#L34-L44

process PAVE_GERMLINE {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/pave:1.5--0'

    input:
    tuple val(meta), path(sage_vcf)
    path genome_fasta
    val genome_ver
    path genome_fai
    path sage_blocklist_regions
    path sage_blocklist_sites
    path clinvar_annotations
    path segment_mappability
    path driver_gene_panel
    path ensembl_data_resources
    path gnomad_resource

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: index
    path 'versions.yml'                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def gnomad_args
    if (genome_ver == '37') {
        gnomad_args = "-gnomad_freq_file ${gnomad_resource}"
    } else if (genome_ver == '38') {
        gnomad_args = "-gnomad_freq_dir ${gnomad_resource} -gnomad_load_chr_on_demand"
    } else {
        log.error "got bad genome version: ${genome_ver}"
        System.exit(1)
    }

    """
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -sample ${meta.sample_id} \\
            -vcf_file ${sage_vcf} \\
            -ref_genome ${genome_fasta} \\
            -ref_genome_version ${genome_ver} \\
            -clinvar_vcf ${clinvar_annotations} \\
            -driver_gene_panel ${driver_gene_panel} \\
            -mappability_bed ${segment_mappability} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -blacklist_bed ${sage_blocklist_regions} \\
            -blacklist_vcf ${sage_blocklist_sites} \\
            -gnomad_pon_filter -1 \\
            ${gnomad_args} \\
            -read_pass_only \\
            -output_dir ./

    # NOTE(SW): hard coded since there is no reliable way to obtain version information.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pave: 1.5
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.sage.pave_germline.vcf.gz{,.tbi}
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
