process NEO_FINDER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-neo:1.2_beta--hdfd78af_0' :
        'biocontainers/hmftools-neo:1.2_beta--hdfd78af_0' }"

    input:
    tuple val(meta), path(purple_dir), path(linx_annotation_dir)
    path genome_fasta
    val genome_ver
    path genome_fai
    path ensembl_data_resources

    output:
    tuple val(meta), path('neo_finder/'), emit: neo_finder_dir
    path 'versions.yml'                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p neo_finder/

    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -sample ${meta.sample_id} \\
            -linx_dir ${linx_annotation_dir} \\
            -somatic_vcf ${purple_dir}/${meta.sample_id}.purple.somatic.vcf.gz \\
            -ref_genome ${genome_fasta} \\
            -ref_genome_version ${genome_ver} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -log_debug \\
            -output_dir neo_finder/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        neo: \$(java -jar ${task.ext.jarPath} -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p neo_finder/
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}

