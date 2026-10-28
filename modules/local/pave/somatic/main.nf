process PAVE_SOMATIC {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-pave:1.9--hdfd78af_0' :
        'biocontainers/hmftools-pave:1.9--hdfd78af_0' }"

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
    val sequencing_platform

    output:
    tuple val(meta), path('pave_somatic/')                  , topic: pave_somatic_dir
    tuple val(meta), val('pave_somatic'), path('.command.*'), topic: command_files
    path 'versions.yml'                                     , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def gnomad_args
    if (genome_ver == '37') {
        gnomad_args = "-gnomad_freq_file ${gnomad_resource}"
    } else if (genome_ver == '38') {
        gnomad_args = "-gnomad_freq_dir ${gnomad_resource}"
    } else {
        error "got bad genome version: ${genome_ver}"
    }

    // Targeted mode
    def pon_artefact_arg = pon_artefacts ? "-pon_artefact_file ${pon_artefacts}" : ''

    """
    mkdir -p pave_somatic/

    pave \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -input_vcf ${sage_vcf} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        ${pon_artefact_arg} \\
        -pon_file ${sage_pon} \\
        ${gnomad_args} \\
        -clinvar_vcf ${clinvar_annotations} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -mappability_bed ${segment_mappability} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -sequencing_type ${sequencing_platform.toUpperCase()} \\
        -threads ${task.cpus} \\
        ${log_level_arg} \\
        -output_dir pave_somatic/ \\
        -output_vcf pave_somatic/${meta.sample_id}.pave.somatic.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pave: \$(pave -version | sed -n '/^Pave version / { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p pave_somatic/

    gzip <<< '' > pave_somatic/${meta.sample_id}.pave.somatic.vcf.gz
    touch pave_somatic/${meta.sample_id}.pave.somatic.vcf.gz.tbi

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
