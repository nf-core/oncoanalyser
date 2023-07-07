// NOTE(SW): use of tumor sample name here is consistent with Pipeline5
//  - https://github.com/hartwigmedical/pipeline5/blob/v5.32/cluster/src/main/java/com/hartwig/pipeline/tertiary/pave/PaveGermline.java#L35-L39
//  - https://github.com/hartwigmedical/pipeline5/blob/v5.32/cluster/src/main/java/com/hartwig/pipeline/tertiary/pave/PaveArguments.java#L34-L44

process PAVE_GERMLINE {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/pave:1.4.3--0'

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
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -sample ${meta.id} \\
            -vcf_file ${sage_vcf} \\
            -ref_genome ${genome_fasta} \\
            -ref_genome_version ${genome_ver} \\
            -clinvar_vcf ${clinvar_annotations} \\
            -driver_gene_panel ${driver_gene_panel} \\
            -mappability_bed ${segment_mappability} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -blacklist_bed ${sage_blocklist_regions} \\
            -blacklist_vcf ${sage_blocklist_sites} \\
            -read_pass_only \\
            -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pave: \$(java -jar ${task.ext.jarPath} 2>&1 | sed -n 's/^.*version: //p')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sage.pave_germline.vcf.gz{,.tbi}
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
