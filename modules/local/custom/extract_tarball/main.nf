process CUSTOM_EXTRACTTARBALL {
    label 'process_single'

    conda "conda-forge::tar=1.34"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'quay.io/nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(tarball)

    output:
    path "${meta.id}/", emit: extracted_dir

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p ${meta.id}/

    tar ${args} -xzvf ${tarball} --strip-components 1 -C ${meta.id}/
    """

    stub:
    """
    mkdir -p ${meta.id}/
    """
}
