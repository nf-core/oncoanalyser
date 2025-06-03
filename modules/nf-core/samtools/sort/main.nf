process SAMTOOLS_SORT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"),  emit: bam
    path 'command.*.{sh,out,err}' ,  emit: logs
    path "versions.yml"           ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def prefix = task.ext.prefix ?: "${meta.prefix}"
    if ("${bam}" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    def log_file_id = "${task.process.split(':')[-1]}.${prefix}"

    """
    samtools sort \\
        ${args} \\
        -T ${prefix} \\
        --threads ${task.cpus} \\
        -o ${prefix}.bam \\
        ${bam}

    for log_file_ext in sh out err; do
        cp .command.\${log_file_ext} command.${log_file_id}.\${log_file_ext}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.prefix}"

    """
    touch ${prefix}.bam

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
