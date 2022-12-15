// NOTE(SW): PAVE gnomad filtering is not yet documented but is used in Pipeline5:
//  - cluster/src/main/java/com/hartwig/pipeline/tertiary/pave/PaveArguments.java#L27-L28
// NOTE(SW): use of tumor sample name here is consistent with pipeline5
//  - cluster/src/main/java/com/hartwig/pipeline/tertiary/pave/PaveGermline.java#L35-L39
//  - cluster/src/main/java/com/hartwig/pipeline/tertiary/pave/PaveArguments.java#L34-L44

process PAVE_GERMLINE {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/pave:1.2.2--0'

    input:
    tuple val(meta), path(sage_vcf)
    path genome_fasta
    path genome_fai
    val genome_ver
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
    java \\
        -Xmx${task.memory.giga}g \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -sample ${meta.id} \\
            -ref_genome_version ${genome_ver} \\
            -ref_genome ${genome_fasta} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -driver_gene_panel ${driver_gene_panel} \\
            -clinvar_annotations ${clinvar_annotations} \\
            -blacklist_bed ${sage_blocklist_regions} \\
            -blacklist_vcf ${sage_blocklist_sites} \\
            -segment_mappability ${segment_mappability} \\
            -vcf_file ${sage_vcf} \\
            -read_pass_only \\
            -output_dir ./

    # NOTE(SW): hard coded since there is no reliable way to obtain version information.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pave: 1.2.2
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sage.pave_germline.vcf.gz{,.tbi}
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
