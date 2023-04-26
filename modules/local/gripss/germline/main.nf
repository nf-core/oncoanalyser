process GRIPSS_GERMLINE {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/gripss:2.3.2--0'

    input:
    tuple val(meta), path(gridss_vcf)
    path genome_fasta
    path genome_fai
    val genome_ver
    path pon_breakends
    path pon_breakpoints
    path known_fusions
    path repeatmasker_annotations

    output:
    tuple val(meta), path('*gripss.filtered.germline.vcf.gz'), path('*gripss.filtered.germline.vcf.gz.tbi'), emit: vcf
    tuple val(meta), path('*.gripss.germline.vcf.gz'), path('*.gripss.germline.vcf.gz.tbi')                , emit: vcf_unfiltered
    path 'versions.yml'                                                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -sample ${meta.normal_id} \\
            -ref_genome_version ${genome_ver} \\
            -ref_genome ${genome_fasta} \\
            -pon_sgl_file ${pon_breakends} \\
            -pon_sv_file ${pon_breakpoints} \\
            -known_hotspot_file ${known_fusions} \\
            -repeat_mask_file ${repeatmasker_annotations} \\
            -vcf ${gridss_vcf} \\
            -output_id germline \\
            -output_dir ./

    # NOTE(SW): hard coded since there is no reliable way to obtain version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gripss: 2.3.2
    END_VERSIONS
    """

    stub:
    """
    cat <<EOF > ${meta.normal_id}.gripss.filtered.germline.vcf.gz
    ##fileformat=VCFv4.1
    ##contig=<ID=.>
    #CHROM      POS     ID      REF     ALT     QUAL    FILTER  INFO
    .   .       .       .       .       .       .
    EOF
    touch ${meta.normal_id}.gripss.filtered.germline.vcf.gz.tbi
    touch ${meta.normal_id}.gripss.germline.vcf.gz
    touch ${meta.normal_id}.gripss.germline.vcf.gz.tbi
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
