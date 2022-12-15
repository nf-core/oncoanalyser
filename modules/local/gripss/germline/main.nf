process GRIPSS_GERMLINE {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/gripss:2.1--0'

    input:
    tuple val(meta), path(gridss_vcf)
    path genome_fasta
    path genome_fai
    val genome_ver
    path pon_breakends
    path pon_breakpoints
    path known_fusions

    output:
    tuple val(meta), path('*.gripss.filtered.vcf.gz'), path('*.gripss.filtered.vcf.gz.tbi'), emit: vcf_hard
    tuple val(meta), path('*.gripss.vcf.gz'), path('*.gripss.vcf.gz.tbi')                  , emit: vcf_soft
    path 'versions.yml'                                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    java \\
        -Xmx${task.memory.giga}g \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -sample ${meta.normal_id} \\
            -ref_genome_version ${genome_ver} \\
            -ref_genome ${genome_fasta} \\
            -pon_sgl_file ${pon_breakends} \\
            -pon_sv_file ${pon_breakpoints} \\
            -known_hotspot_file ${known_fusions} \\
            -vcf ${gridss_vcf} \\
            -output_dir ./

    # NOTE(SW): hard coded since there is no reliable way to obtain version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gripss: 2.1
    END_VERSIONS
    """

    stub:
    """
    cat <<EOF > ${meta.normal_id}.gripss.filtered.vcf.gz
    ##fileformat=VCFv4.1
    ##contig=<ID=.>
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    .	.	.	.	.	.	.
    EOF
    touch ${meta.normal_id}.gripss.filtered.vcf.gz.tbi
    touch ${meta.normal_id}.gripss.vcf.gz
    touch ${meta.normal_id}.gripss.vcf.gz.tbi
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
