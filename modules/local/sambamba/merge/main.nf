process SAMBAMBA_MERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sambamba:1.0.1--h6f6fda4_0' :
        'biocontainers/sambamba:1.0.1--h6f6fda4_0' }"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path('*bam'), emit: bam
    path 'command.*.{sh,out,err}', emit: logs
    path 'versions.yml'          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def log_file_id = "${task.process.split(':')[-1]}.${meta.sample_id}"

    """
    sambamba merge \\
        ${args} \\
        --nthreads ${task.cpus} \\
        ${meta.sample_id}.bam \\
        ${bams}

    for log_file_ext in sh out err; do
        cp .command.\${log_file_ext} command.${log_file_id}.\${log_file_ext}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sambamba: \$(sambamba --version 2>&1 | sed -n '/^sambamba / { s/^.* //p }' | head -n1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}.bam

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
