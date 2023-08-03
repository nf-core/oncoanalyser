process GPGR_LINX {
    tag "${meta.id}"
    label 'process_single'

    container 'docker.io/scwatts/gpgr:1.4.5'

    input:
    tuple val(meta), path(linx_annotation_dir), path(linx_visualiser_dir)

    output:
    tuple val(meta), path('*_linx.html'), emit: html
    path 'versions.yml'                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    gpgr.R linx \\
        ${args} \\
        --sample ${meta.id} \\
        --plot ${linx_visualiser_dir}/ \\
        --table ${linx_annotation_dir}/ \\
        --out ${meta.id}_linx.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/^R version \\([0-9.]\\+\\).\\+/\\1/')
        gpgr: \$(gpgr.R --version | cut -f2 -d' ')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_linx.html
    echo -e '${task.process}:\n  stub: noversions\n' > versions.yml
    """
}
