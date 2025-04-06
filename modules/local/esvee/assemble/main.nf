process ESVEE_ASSEMBLE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-esvee:1.0.3--hdfd78af_0' :
        'biocontainers/hmftools-esvee:1.0.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_prep_bam), path(tumor_prep_bai), path(normal_prep_bam), path(normal_prep_bai), path(prep_dir)
    path genome_fasta
    path genome_fai
    path genome_dict
    path genome_img
    val genome_ver
    path decoy_sequences_image

    output:
    tuple val(meta), path('assemble/')                                 , emit: assemble_dir
    tuple val(meta), path("assemble/${meta.tumor_id}.esvee.raw.vcf.gz"), emit: raw_vcf
    path 'versions.yml'                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def reference_arg = meta.normal_id != null ? "-reference ${meta.normal_id}" : ''
    def reference_bam_arg = meta.normal_id != null ? "-reference_bam ${normal_prep_bam}" : ''

    def decoy_genome_arg = decoy_sequences_image ? "-decoy_genome ${decoy_sequences_image}" : ''

    """
    mkdir -p assemble/

    esvee com.hartwig.hmftools.esvee.assembly.AssemblyApplication \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -tumor ${meta.tumor_id} \\
        -tumor_bam ${tumor_prep_bam} \\
        ${reference_arg} \\
        ${reference_bam_arg} \\
        -esvee_prep_dir ${prep_dir}/ \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        ${decoy_genome_arg} \\
        -write_types 'JUNC_ASSEMBLY;PHASED_ASSEMBLY;ALIGNMENT;BREAKEND;VCF' \\
        -output_dir assemble/ \\
        -threads ${task.cpus} \\
        -perf_log_time 10 \\
        -log_level DEBUG

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esvee: \$(esvee -version | sed -n '/^Esvee version/ { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p assemble/

    touch assemble/${meta.tumor_id}.esvee.raw.vcf.gz
    touch assemble/${meta.tumor_id}.esvee.raw.vcf.gz.tbi
    touch assemble/${meta.tumor_id}.esvee.alignment.tsv
    touch assemble/${meta.tumor_id}.esvee.assembly.tsv
    touch assemble/${meta.tumor_id}.esvee.phased_assembly.tsv
    touch assemble/${meta.tumor_id}.esvee.breakend.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
