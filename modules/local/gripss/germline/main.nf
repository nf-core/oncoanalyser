process GRIPSS_GERMLINE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-gripss:2.4--hdfd78af_0' :
        'quay.io/biocontainers/hmftools-gripss:2.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(gridss_vcf)
    path genome_fasta
    val genome_ver
    path genome_fai
    path pon_breakends
    path pon_breakpoints
    path known_fusions
    path repeatmasker_annotations

    output:
    tuple val(meta), path('*.filtered.germline.vcf.gz'), path('*.filtered.germline.vcf.gz.tbi'), emit: vcf
    tuple val(meta), path('*gripss.germline.vcf.gz'), path('*gripss.germline.vcf.gz.tbi')      , emit: vcf_unfiltered
    path 'versions.yml'                                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    gripss \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        -sample ${meta.normal_id} \\
        -reference ${meta.tumor_id} \\
        -vcf ${gridss_vcf} \\
        -germline \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -pon_sgl_file ${pon_breakends} \\
        -pon_sv_file ${pon_breakpoints} \\
        -known_hotspot_file ${known_fusions} \\
        -repeat_mask_file ${repeatmasker_annotations} \\
        -output_id germline \\
        -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gripss: \$(gripss -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.normal_id}.gripss.filtered.germline.vcf.gz
    touch ${meta.normal_id}.gripss.filtered.germline.vcf.gz.tbi
    touch ${meta.normal_id}.gripss.germline.vcf.gz
    touch ${meta.normal_id}.gripss.germline.vcf.gz.tbi
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
