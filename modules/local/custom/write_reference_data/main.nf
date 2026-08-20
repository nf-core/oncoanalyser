nextflow.enable.types = true

process WRITE_REFERENCE_DATA {
    tag "${fp.name}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:24.04' :
        'nf-core/ubuntu:24.04' }"

    input:
    fp: Path

    topic:
    file(fp, includeInputs: true) >> 'write_reference_data'

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    """
}
