process GRIPSS_SOMATIC {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/gripss:2.3.5--0'

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
    tuple val(meta), path('*.filtered.somatic.vcf.gz'), path('*.filtered.somatic.vcf.gz.tbi'), emit: vcf
    tuple val(meta), path('*.gripss.somatic.vcf.gz'), path('*.gripss.somatic.vcf.gz.tbi')    , emit: vcf_unfiltered
    path 'versions.yml'                                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -sample ${meta.tumor_id} \\
            -reference ${meta.normal_id} \\
            -vcf ${gridss_vcf} \\
            -ref_genome ${genome_fasta} \\
            -ref_genome_version ${genome_ver} \\
            -pon_sgl_file ${pon_breakends} \\
            -pon_sv_file ${pon_breakpoints} \\
            -known_hotspot_file ${known_fusions} \\
            -repeat_mask_file ${repeatmasker_annotations} \\
            -output_id somatic \\
            -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gripss: \$(java -jar ${task.ext.jarPath} | sed -n '1s/^.*version: //p')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.tumor_id}.gripss.filtered.somatic.vcf.gz
    touch ${meta.tumor_id}.gripss.filtered.somatic.vcf.gz.tbi
    touch ${meta.tumor_id}.gripss.somatic.vcf.gz
    touch ${meta.tumor_id}.gripss.somatic.vcf.gz.tbi
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
