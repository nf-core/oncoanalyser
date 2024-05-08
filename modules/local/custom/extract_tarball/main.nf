process CUSTOM_EXTRACTTARBALL {
    label 'process_single'

    conda "conda-forge::tar"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'quay.io/nf-core/ubuntu:20.04' }"

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
