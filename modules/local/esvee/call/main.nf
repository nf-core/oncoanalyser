process ESVEE_CALL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sv-prep:1.2.4--hdfd78af_0' :
        'docker.io/scwatts/hmftools-esvee:1.0_beta--hdfd78af_0--1' }"

    input:
    tuple val(meta), path(ref_depth_vcf), path(fragment_lengths_tsv)
    path genome_fasta
    val genome_ver
    path pon_breakends
    path pon_breakpoints
    path known_fusions
    path repeatmasker_annotations

    output:
    tuple val(meta), path("caller/"), emit: caller_dir
    tuple val(meta), path("caller/${meta.tumor_id}.esvee.germline.vcf.gz"),   path("caller/${meta.tumor_id}.esvee.germline.vcf.gz.tbi"),   emit: germline_vcf
    tuple val(meta), path("caller/${meta.tumor_id}.esvee.somatic.vcf.gz"),    path("caller/${meta.tumor_id}.esvee.somatic.vcf.gz.tbi"),    emit: somatic_vcf
    tuple val(meta), path("caller/${meta.tumor_id}.esvee.unfiltered.vcf.gz"), path("caller/${meta.tumor_id}.esvee.unfiltered.vcf.gz.tbi"), emit: unfiltered_vcf
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p caller/

    # Esvee expects the fragment_lengths.tsv input file to be in `output_dir`
    ln -sf \$(realpath ${fragment_lengths_tsv}) caller/

    esvee com.hartwig.hmftools.esvee.caller.CallerApplication \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        -sample ${meta.tumor_id} \\
        -reference ${meta.normal_id} \\
        -input_vcf ${ref_depth_vcf} \\
        -ref_genome_version ${genome_ver} \\
        -known_hotspot_file ${known_fusions} \\
        -pon_sgl_file ${pon_breakends} \\
        -pon_sv_file ${pon_breakpoints} \\
        -repeat_mask_file ${repeatmasker_annotations} \\
        -output_dir caller/ \\
        -log_debug

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esvee: \$(esvee -version | sed 's/^.* //')
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

    echo \${vcf_template} | gzip -c > caller/${meta.tumor_id}.esvee.somatic.vcf.gz
    echo \${vcf_template} | gzip -c > caller/${meta.tumor_id}.esvee.germline.vcf.gz
    echo \${vcf_template} | gzip -c > caller/${meta.tumor_id}.esvee.unfiltered.vcf.gz

    touch caller/${meta.tumor_id}.esvee.somatic.vcf.gz.tbi
    touch caller/${meta.tumor_id}.esvee.germline.vcf.gz.tbi
    touch caller/${meta.tumor_id}.esvee.unfiltered.vcf.gz.tbi

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}