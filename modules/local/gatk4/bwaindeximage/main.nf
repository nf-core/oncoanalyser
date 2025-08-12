process GATK4_BWA_INDEX_IMAGE {
    tag "${genome_fasta.name}"
    label 'process_medium'

    conda "bioconda::gatk4=4.6.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.6.1.0--py310hdfd78af_0' :
        'biocontainers/gatk4:4.6.1.0--py310hdfd78af_0' }"

    input:
    path genome_fasta

    output:
    path "${genome_fasta}.img", emit: img
    path 'versions.yml'       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    gatk \\
        BwaMemIndexImageCreator \\
        ${args} \\
        -I ${genome_fasta} \\
        -O ${genome_fasta}.img

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version | sed -n '/GATK/ { s/^.* v//p }')
    END_VERSIONS
    """

    stub:
    """
    touch ${genome_fasta}.img

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
