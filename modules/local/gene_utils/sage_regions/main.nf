process GENE_UTILS_SAGE_REGIONS {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-gene-utils:1.2--hdfd78af_0' :
        'biocontainers/hmftools-gene-utils:1.2--hdfd78af_0' }"

    input:
    path driver_gene_panel
    val genome_ver
    path clinvar_annotations
    path ensembl_data_resources

    output:
    path 'ActionableCodingPanel.*.bed.gz'
    path 'CoverageCodingPanel.*.bed.gz'
    path 'versions.yml'                  , emit: versions
    path '.command.{sh,log}'             , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    """
    mkdir -p resources/{ensembl_data_cache,sage}/${genome_ver}/

    ln -s \$(pwd)/${clinvar_annotations} resources/sage/${genome_ver}/
    ln -s \$(pwd)/${ensembl_data_resources}/* resources/ensembl_data_cache/${genome_ver}/

    gene-utils \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        com.hartwig.hmftools.geneutils.drivers.GenerateDriverGeneFiles \\
        -ref_genome_version ${genome_ver} \\
        -resource_repo_dir resources/ \\
        -driver_gene_panel ${driver_gene_panel} \\
        -log_debug \\
        -output_dir ./ \\
        -log_level ${params.module_log_level}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gene-utils: \$(gene-utils com.hartwig.hmftools.geneutils.drivers.GenerateDriverGeneFiles -version | sed -n '/GeneUtils/ { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    touch ActionableCodingPanel.${genome_ver}.bed.gz
    touch CoverageCodingPanel.${genome_ver}.bed.gz

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
