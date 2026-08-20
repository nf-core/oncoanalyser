// NOTE(SW): kept untyped (legacy `topic:` qualifier) because this process is aliased multiple times
// in prepare_reference, and a typed process with a `topic:` section under multiple aliases hangs the
// run (Nextflow issue #7434). Re-enable static types once that bug is fixed.
process WRITE_REFERENCE_DATA {
    tag "${fp.name}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:24.04' :
        'nf-core/ubuntu:24.04' }"

    input:
    path fp

    output:
    path fp, includeInputs: true, topic: write_reference_data

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    """
}
