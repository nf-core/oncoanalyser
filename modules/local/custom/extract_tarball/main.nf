nextflow.enable.types = true

process CUSTOM_EXTRACTTARBALL {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:24.04' :
        'nf-core/ubuntu:24.04' }"

    input:
    tuple(meta: Map, tarball: Path)

    topic:
    tuple(meta, file("${meta.id}/")) >> 'extracted_dir'
    tuple(meta, 'extracttarball', files('.command.*')) >> 'command_files'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p ${meta.id}/

    tar ${args} -xzvf ${tarball} --transform 's#^\\./##' --strip-components 1 -C ${meta.id}/
    """

    stub:
    """
    mkdir -p ${meta.id}/
    """
}
