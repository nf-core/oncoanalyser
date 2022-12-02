// NOTE(SW): PAVE gnomad filtering is not yet documented but is used in Pipeline5:
//  - cluster/src/main/java/com/hartwig/pipeline/tertiary/pave/PaveArguments.java#L27-L28
// NOTE(SW): use of tumor sample name here is consistent with pipeline5
//  - cluster/src/main/java/com/hartwig/pipeline/tertiary/pave/PaveGermline.java#L35-L39
//  - cluster/src/main/java/com/hartwig/pipeline/tertiary/pave/PaveArguments.java#L34-L44

process PAVE_GERMLINE {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/pave:1.4--0'

    input:
    tuple val(meta), path(sage_vcf)
    path genome_fasta
    path genome_fai
    val genome_ver
    path sage_blacklist_bed
    path sage_blacklist_vcf
    path clinvar_vcf
    path mappability_bed
    path driver_gene_panel
    path ensembl_data_dir

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
            -ensembl_data_dir ${ensembl_data_dir} \\
            -driver_gene_panel ${driver_gene_panel} \\
            -clinvar_vcf ${clinvar_vcf} \\
            -blacklist_bed ${sage_blacklist_bed} \\
            -blacklist_vcf ${sage_blacklist_vcf} \\
            -mappability_bed ${mappability_bed} \\
            -vcf_file ${sage_vcf} \\
            -read_pass_only \\
            -output_dir ./

    # NOTE(SW): hard coded since there is no reliable way to obtain version information.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pave: 1.4
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sage.pave_germline.vcf.gz{,.tbi}
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
