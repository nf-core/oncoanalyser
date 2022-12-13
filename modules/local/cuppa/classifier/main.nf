process CUPPA_CLASSIFIER {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/cuppa:1.6--0'

    input:
    tuple val(meta), path(isofox_dir), path(purple_dir), path(linx_dir), path(virusinterpreter_dir)
    path reference_data

    output:
    tuple val(meta), path('*csv'), emit: csv
    path 'versions.yml'          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def categories_arg = isofox_dir ? 'ALL' : 'DNA'

    """
    # Symlink input files into a single directory
    mkdir sample_data/
    input_dirs="${isofox_dir} ${purple_dir} ${linx_dir} ${virusinterpreter_dir}"
    find -L \${input_dirs} -maxdepth 1 -type f -exec ln -s ../{} sample_data/ \\;

    java \\
        -Xmx${task.memory.giga}g \\
        -jar ${task.ext.jarPath} \\
            -categories ${categories_arg} \\
            -ref_data_dir ${reference_data} \\
            -sample_data ${meta.id} \\
            -sample_data_dir sample_data/ \\
            -output_dir ./

    # NOTE(SW): hard coded since there is no reliable way to obtain version information.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cuppa: 1.6
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.cup.data.csv
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
