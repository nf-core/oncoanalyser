process GRIPSS_SOMATIC {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/gripss:2.3.1--0'

    input:
    tuple val(meta), path(gridss_vcf)
    path genome_fasta
    path genome_fai
    val genome_ver
    path breakend_pon
    path breakpoint_pon
    path known_fusions
    path repeat_mask_file

    output:
    tuple val(meta), path('*.gripss.filtered.vcf.gz'), path('*.gripss.filtered.vcf.gz.tbi'), emit: vcf_hard
    tuple val(meta), path('*.gripss.vcf.gz'), path('*.gripss.vcf.gz.tbi')                  , emit: vcf_soft
    path 'versions.yml'                                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def repeat_mask_file_arg = repeat_mask_file ? "-repeat_mask_file ${repeat_mask_file}" : ''

    """
    java \\
        -Xmx${task.memory.giga}g \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -sample ${meta.tumor_id} \\
            -reference ${meta.normal_id} \\
            -ref_genome_version ${genome_ver} \\
            -ref_genome ${genome_fasta} \\
            -pon_sgl_file ${breakend_pon} \\
            -pon_sv_file ${breakpoint_pon} \\
            -known_hotspot_file ${known_fusions} \\
            -vcf ${gridss_vcf} \\
            ${repeat_mask_file_arg} \\
            -output_dir ./

    # NOTE(SW): hard coded since there is no reliable way to obtain version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gripss: 2.3.1
    END_VERSIONS
    """

    stub:
    """
    cat <<EOF > ${meta.tumor_id}.gripss.filtered.vcf.gz
    ##fileformat=VCFv4.1
    ##contig=<ID=.>
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    .	.	.	.	.	.	.
    EOF
    touch ${meta.tumor_id}.gripss.filtered.vcf.gz.tbi
    touch ${meta.tumor_id}.gripss.vcf.gz
    touch ${meta.tumor_id}.gripss.vcf.gz.tbi
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
