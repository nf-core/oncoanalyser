process CUSTOM_EXTRACTTARBALL {
    label 'process_single'

    conda "conda-forge::tar"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/debian:bullseye-slim' :
        'docker pull debian:bullseye-slim' }"

    input:
    tuple val(meta), path(tarball)

    output:
    path "${meta.id}/", emit: dir

    script:
    """
    mkdir -p ${meta.id}/
    tar -xzvf ${tarball} --strip-components 1 -C ${meta.id}/
    """

    stub:
    """
    mkdir -p ${meta.id}/
    """
}
