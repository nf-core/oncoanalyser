process GRIPSS_SOMATIC {
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
    path target_region_bed

    output:
    tuple val(meta), path('*.gripss.filtered{,.somatic}.vcf.gz'), path('*.gripss.filtered{,.somatic}.vcf.gz.tbi'), emit: vcf
    tuple val(meta), path('*.gripss{,.somatic}.vcf.gz'), path('*.gripss{,.somatic}.vcf.gz.tbi')                  , emit: vcf_unfiltered
    path 'versions.yml'                                                                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def reference_arg = meta.containsKey('normal_id') ? "-reference ${meta.normal_id}" : ''
    def target_regions_bed_arg = target_region_bed ? "-target_regions_bed ${target_region_bed}" : ''
    def output_id_arg = meta.containsKey('normal_id') ? '-output_id somatic' : ''

    """
    gripss \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        -sample ${meta.tumor_id} \\
        ${reference_arg} \\
        -vcf ${gridss_vcf} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -pon_sgl_file ${pon_breakends} \\
        -pon_sv_file ${pon_breakpoints} \\
        -known_hotspot_file ${known_fusions} \\
        -repeat_mask_file ${repeatmasker_annotations} \\
        ${target_regions_bed_arg} \\
        ${output_id_arg} \\
        -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gripss: \$(gripss -version | sed 's/^.* //')
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
