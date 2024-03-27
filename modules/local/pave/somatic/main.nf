process PAVE_SOMATIC {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-pave:1.6--hdfd78af_0' :
        'quay.io/biocontainers/hmftools-pave:1.6--hdfd78af_0' }"

    input:
    tuple val(meta), path(sage_vcf), path(sage_tbi)
    path genome_fasta
    val genome_ver
    path genome_fai
    path sage_pon
    path pon_artefacts
    path sage_blocklist_regions
    path sage_blocklist_sites
    path clinvar_annotations
    path segment_mappability
    path driver_gene_panel
    path ensembl_data_resources
    path gnomad_resource

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: index
    path 'versions.yml'                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def pon_filters
    def gnomad_args
    if (genome_ver == '37') {
        pon_filters = 'HOTSPOT:10:5;PANEL:6:5;UNKNOWN:6:0'
        gnomad_args = "-gnomad_freq_file ${gnomad_resource}"
    } else if (genome_ver == '38') {
        pon_filters = 'HOTSPOT:5:5;PANEL:2:5;UNKNOWN:2:0'
        gnomad_args = "-gnomad_freq_dir ${gnomad_resource}"
    } else {
        log.error "got bad genome version: ${genome_ver}"
        System.exit(1)
    }

    // Targeted mode
    def pon_artefact_arg = pon_artefacts ? "-pon_artefact_file ${pon_artefacts}" : ''
    def sage_blocklist_regions_arg = sage_blocklist_regions ? "-blacklist_bed ${sage_blocklist_regions}" : ''
    def sage_blocklist_sites_arg = sage_blocklist_sites ? "-blacklist_vcf ${sage_blocklist_sites}" : ''
    def clinvar_annotations = clinvar_annotations ? "-clinvar_vcf ${clinvar_annotations}" : ''

    """
    pave \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -sample ${meta.sample_id} \\
        -vcf_file ${sage_vcf} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_ver} \\
        -pon_file ${sage_pon} \\
        -pon_filters "${pon_filters}" \\
        ${pon_artefact_arg} \\
        ${clinvar_annotations} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -mappability_bed ${segment_mappability} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        ${sage_blocklist_regions_arg} \\
        ${sage_blocklist_sites_arg} \\
        ${gnomad_args} \\
        -read_pass_only \\
        -threads ${task.cpus} \\
        -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pave: \$(pave -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.sage.pave_somatic.vcf.gz{,.tbi}
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
