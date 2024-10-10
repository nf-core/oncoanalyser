process ESVEE_ASSEMBLE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-esvee:1.0_beta--hdfd78af_0' :
        'biocontainers/hmftools-esvee:1.0_beta--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_prep_bam), path(normal_prep_bam), path(junctions_tsv), path(fragment_lengths_tsv)
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

    def reference_arg = meta.normal_id != null ? "-reference ${meta.normal_id}" : ""
    def reference_bam_arg = meta.normal_id != null ? "-reference_bam ${normal_prep_bam}" : ""

    def decoy_genome_arg = decoy_sequences_image ? "-decoy_genome ${decoy_sequences_image}" : ""

    """
    mkdir -p assemble/

    # Esvee expects the fragment_lengths.tsv input file to be in `output_dir`
    ln -sf \$(realpath ${fragment_lengths_tsv}) assemble/

    esvee com.hartwig.hmftools.esvee.EsveeApplication \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        -tumor ${meta.tumor_id} \\
        -tumor_bam ${tumor_prep_bam} \\
        ${reference_arg} \\
        ${reference_bam_arg} \\
        -junction_files ${junctions_tsv} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        ${decoy_genome_arg} \\
        -write_types "JUNC_ASSEMBLY;PHASED_ASSEMBLY;ALIGNMENT;BREAKEND;VCF" \\
        -output_dir assemble/ \\
        -threads ${task.cpus} \\
        -perf_log_time 10 \\
        -log_debug

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esvee: \$(esvee -version | sed 's/^.* //')
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
