process PAVE_SOMATIC {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-pave:1.8.2--hdfd78af_0' :
        'biocontainers/hmftools-pave:1.8.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(sage_vcf), path(sage_tbi)
    path genome_fasta
    val genome_ver
    path genome_fai
    path pon_artefacts
    path sage_pon
    path clinvar_annotations
    path segment_mappability
    path driver_gene_panel
    path ensembl_data_resources
    path gnomad_resource
    val sequencing_type

    output:
    tuple val(meta), path('*.vcf.gz')    , emit: vcf
    tuple val(meta), path('*.vcf.gz.tbi'), emit: index
    path 'versions.yml'                  , emit: versions
    path '.command.*'                    , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def gnomad_arg = genome_ver == '38'
        ? "-gnomad_freq_dir ${gnomad_resource}"
        : "-gnomad_freq_file ${gnomad_resource}"

    // Targeted mode
    def pon_artefact_arg = pon_artefacts ? "-pon_artefact_file ${pon_artefacts}" : ''

    """
    pave \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -input_vcf ${sage_vcf} \\
        -output_vcf ${meta.sample_id}.pave.somatic.vcf.gz \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        ${pon_artefact_arg} \\
        -pon_file ${sage_pon} \\
        ${gnomad_arg} \\
        -clinvar_vcf ${clinvar_annotations} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -mappability_bed ${segment_mappability} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -sequencing_type ${sequencing_type} \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pave: \$(pave -version | sed -n '/^Pave version / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.pave.somatic.vcf.gz{,.tbi}

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
