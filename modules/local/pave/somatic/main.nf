// NOTE(SW): PAVE gnomad filtering is not yet documented but is used in Pipeline5 https://github.com/hartwigmedical/pipeline5/blob/master/cluster/src/main/java/com/hartwig/pipeline/tertiary/pave/PaveArguments.java#L27-L28

process PAVE_SOMATIC {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/pave:1.2.2--0'

    input:
    tuple val(meta), path(sage_vcf)
    path genome_fasta
    path genome_fai
    val genome_ver
    path sage_pon
    path segment_mappability
    path driver_gene_panel
    path ensembl_data_resources
    path gnomad_pon_dir

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: index
    path 'versions.yml'                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def pon_filters
    def gnomad_args
    if (genome_ver == '37') {
        pon_filters = 'HOTSPOT:10:5;PANEL:6:5;UNKNOWN:6:0'
        gnomad_args = ''
    } else if (genome_ver == '38') {
        pon_filters = 'HOTSPOT:5:5;PANEL:2:5;UNKNOWN:2:0'
        gnomad_args = "-gnomad_freq_dir \"${gnomad_pon_dir}\" -gnomad_load_chr_on_demand"
    } else {
        log.error "got bad genome version: ${genome_ver}"
        System.exit(1)
    }

    """
    java \\
        -Xmx${task.memory.giga}g \\
        -jar ${task.ext.jarPath} \\
            -sample ${meta.id} \\
            -ref_genome_version ${genome_ver} \\
            -ref_genome ${genome_fasta} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -driver_gene_panel ${driver_gene_panel} \\
            -pon_file ${sage_pon} \\
            -pon_filters "${pon_filters}" \\
            -mappability_bed ${segment_mappability} \\
            -vcf_file ${sage_vcf} \\
            -read_pass_only \\
            -write_pass_only \\
            ${gnomad_args} \\
            -output_dir ./

    # NOTE(SW): hard coded since there is no reliable way to obtain version information.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pave: 1.2.2
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sage.pave_somatic.vcf.gz{,.tbi}
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
