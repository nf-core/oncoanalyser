process CUSTOM_EXTRACTTARBALL {
    label 'process_single'

    container 'docker.io/ubuntu:20.04'

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
