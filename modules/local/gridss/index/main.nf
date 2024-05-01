process GRIDSS_INDEX {
    tag "${genome_fasta.name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h50ea8bc_3' :
        'quay.io/biocontainers/gridss:2.13.2--h50ea8bc_3' }"

    input:
    path genome_fasta
    path genome_fai
    path genome_dict
    path genome_bwa_index

    output:
    path 'gridss_index', emit: index
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    # Symlink BWA indices next to assembly FASTA
    ln -s \$(find -L ${genome_bwa_index} -type f) ./

    # Create indexes
    PrepareReference \\
        ${args} \\
        REFERENCE_SEQUENCE=${genome_fasta} \\
        CREATE_SEQUENCE_DICTIONARY='false' \\
        CREATE_BWA_INDEX_IMAGE='true' \\
        CREATE_GRIDSS_REFERENCE_CACHE='true'

    # Move under single directory for output
    mkdir -p gridss_index/
    mv ${genome_fasta.name}.img gridss_index/
    mv ${genome_fasta.name}.gridsscache gridss_index/

    # Copy across BWA index
    ln -s ../${genome_fasta.name}.{amb,ann,bwt,pac,sa} gridss_index/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: \$(CallVariants --version 2>&1 | sed 's/-gridss\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p gridss_index/
    touch gridss_index/${genome_fasta.name}.{sa,pac,bwt,ann,amb}
    touch gridss_index/${genome_fasta.name}.img
    touch gridss_index/${genome_fasta.name}.gridsscache

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
