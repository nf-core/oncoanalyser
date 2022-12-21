process CUSTOM_EXTRACTTARBALL {
    label 'process_single'

    container 'public.ecr.aws/ubuntu/ubuntu:20.04_stable'

    input:
    tuple val(meta), path(tarball)

    output:
    path "${meta.id}/", emit: dir

    script:
    """
    mkdir -p ${meta.id}/
    tar -xzvf ${tarball} --strip-components 1 -C ${meta.id}/
    """
}
