process ISOFOX_NEO {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/isofox:1.7.1--0'

    input:
    tuple val(meta), path(bam), path(bai)
    val functions
    val read_length
    path genome_fasta
    val genome_ver
    path genome_fai
    path ensembl_data_resources

    output:
    tuple val(meta), path('isofox/'), emit: isofox_dir
    path 'versions.yml'             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p isofox/

    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            ${args} \\
            -sample ${meta.sample_id} \\
            -bam_file ${bam} \\
            -functions NEO_EPITOPES \\
            -neoepitope_file ${neo_finder_dir}/${meta.sample_id}.neo.neo_data.tsv \\
            -read_length ${read_length} \\
            -ref_genome ${genome_fasta} \\
            -ref_genome_version ${genome_ver} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -threads ${task.cpus} \\
            -output_dir isofox/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isofox: \$(java -jar ${task.ext.jarPath} -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p isofox/
    touch isofox/placeholder
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
