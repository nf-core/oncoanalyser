process LINXREPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-linxreport:1.1.0--r44hdfd78af_0' :
        'biocontainers/r-linxreport:1.1.0--r44hdfd78af_0' }"

    input:
    tuple val(meta), path(linx_annotation_dir), path(linx_visualiser_dir)

    output:
    tuple val(meta), path('*_linx.html'), emit: html
    path 'versions.yml'                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def plot_dir = linx_visualiser_dir.resolve('all/').toUriString()

    """
    # Set input plot directory and create it doesn't exist. See the LINX visualiser module for further info.
    if [[ ! -e ${plot_dir} ]]; then
        mkdir -p ${plot_dir};
    fi;

    linxreport.R \\
        ${args} \\
        --sample ${meta.sample_id} \\
        --plot ${plot_dir} \\
        --table ${linx_annotation_dir} \\
        --out ${meta.sample_id}_linx.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -n1 | sed 's/^R version \\([0-9.]\\+\\).\\+/\\1/')
        linxreport: \$(linxreport.R --version)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.sample_id}_linx.html

    echo -e '${task.process}:\n  stub: noversions\n' > versions.yml
    """
}
