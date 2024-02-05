process NEO_FINDER {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/neo:1.1_beta--0'

    input:
    tuple val(meta), path(purple_dir), path(linx_dir)
    path genome_fasta
    val genome_ver
    path ensembl_data_resources

    output:
    tuple val(meta), path('neo/'), emit: neo_finder_dir
    path 'versions.yml'          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -cp ${task.ext.jarPath} \\
        com.hartwig.hmftools.neo.epitope.NeoEpitopeFinder \\
            ${args} \\
            -sample ${meta.sample_id} \\
            -linx_dir ${linx_dir} \\
            -somatic_vcf ${purple_dir}/${meta.sample_id}.purple.somatic.vcf.gz \\
            -ref_genome ${genome_fasta} \\
            -ref_genome_version ${genome_ver} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -output_dir ${output_dir} \\
            -log_debug \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        neo: \$(java -jar ${task.ext.jarPath} -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p neo/
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}

