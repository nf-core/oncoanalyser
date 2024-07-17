process ESVEE_ASSEMBLE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sv-prep:1.2.4--hdfd78af_0' :
        'docker.io/scwatts/hmftools-esvee:1.0_beta--hdfd78af_0--1' }"

    input:
    tuple val(meta), path(tumor_prep_bam), path(normal_prep_bam), path(junctions_tsv), path(fragment_lengths_tsv)
    path genome_fasta
    path genome_fai
    path genome_dict
    val genome_ver
    path decoy_sequences_image

    output:
    tuple val(meta), path("assemble/${meta.tumor_id}.esvee.raw.vcf.gz"), emit: raw_vcf
    tuple val(meta), path('assemble/')                                 , emit: assemble_dir
    path 'versions.yml'                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p assemble/

    # Esvee expects the fragment_lengths.tsv input file to be in `output_dir`
    ln -sf \$(realpath ${fragment_lengths_tsv}) assemble/

    esvee com.hartwig.hmftools.esvee.EsveeApplication \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        -tumor ${meta.tumor_id} \\
        -reference ${meta.normal_id} \\
        -tumor_bam ${tumor_prep_bam} \\
        -reference_bam ${normal_prep_bam} \\
        -junction_files ${junctions_tsv} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -decoy_genome ${decoy_sequences_image} \\
        -write_types "JUNC_ASSEMBLY;ALIGNMENT;ALIGNMENT_DATA;BREAKEND;VCF" \\
        -asm_ref_base_write_max 0 \\
        -phase_process_limit 500 \\
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
    touch assemble/${meta.tumor_id}.esvee.align_detailed.tsv
    touch assemble/${meta.tumor_id}.esvee.alignment.tsv
    touch assemble/${meta.tumor_id}.esvee.assemblies.tsv
    touch assemble/${meta.tumor_id}.esvee.breakend.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
