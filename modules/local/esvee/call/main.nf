process ESVEE_CALL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-esvee:1.0.3--hdfd78af_0' :
        'biocontainers/hmftools-esvee:1.0.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(ref_depth_vcf), path(prep_dir)
    path genome_fasta
    val genome_ver
    path pon_breakends
    path pon_breakpoints
    path known_fusions
    path repeatmasker_annotations

    output:
    tuple val(meta), path("caller/")                                                                                                     , emit: caller_dir
    tuple val(meta), path("caller/${meta.tumor_id}.esvee.unfiltered.vcf.gz"), path("caller/${meta.tumor_id}.esvee.unfiltered.vcf.gz.tbi"), emit: unfiltered_vcf
    tuple val(meta), path("caller/${meta.tumor_id}.esvee.somatic.vcf.gz"),    path("caller/${meta.tumor_id}.esvee.somatic.vcf.gz.tbi")   , emit: somatic_vcf
    tuple val(meta), path("caller/${meta.tumor_id}.esvee.germline.vcf.gz"),   path("caller/${meta.tumor_id}.esvee.germline.vcf.gz.tbi")  , emit: germline_vcf, optional: true
    path 'versions.yml'                                                                                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def reference_arg = meta.normal_id != null ? "-reference ${meta.normal_id}" : ''

    """
    mkdir -p caller/

    esvee com.hartwig.hmftools.esvee.caller.CallerApplication \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${meta.tumor_id} \\
        ${reference_arg} \\
        -input_vcf ${ref_depth_vcf} \\
        -esvee_prep_dir ${prep_dir}/ \\
        -ref_genome_version ${genome_ver} \\
        -known_hotspot_file ${known_fusions} \\
        -pon_sgl_file ${pon_breakends} \\
        -pon_sv_file ${pon_breakpoints} \\
        -repeat_mask_file ${repeatmasker_annotations} \\
        -output_dir caller/ \\
        -log_level DEBUG

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esvee: \$(esvee -version | sed -n '/^Esvee version/ { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p caller/

    vcf_template='##fileformat=VCFv4.1
    ##contig=<ID=.>
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    .	.	.	.	.	.	.
    '

    echo \${vcf_template} | gzip -c > caller/${meta.tumor_id}.esvee.unfiltered.vcf.gz
    echo \${vcf_template} | gzip -c > caller/${meta.tumor_id}.esvee.somatic.vcf.gz

    touch caller/${meta.tumor_id}.esvee.unfiltered.vcf.gz.tbi
    touch caller/${meta.tumor_id}.esvee.somatic.vcf.gz.tbi

    ${ (meta.normal_id != null) ? "touch caller/${meta.tumor_id}.esvee.germline.vcf.gz" : '' }
    ${ (meta.normal_id != null) ? "touch caller/${meta.tumor_id}.esvee.germline.vcf.gz.tbi" : '' }

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
