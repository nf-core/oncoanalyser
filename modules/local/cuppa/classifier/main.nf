process CUPPA_CLASSIFIER {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/scwatts/cuppa:1.8--1'

    input:
    tuple val(meta), path(isofox_dir), path(purple_dir), path(linx_dir), path(virusinterpreter)
    val ref_genome_ver
    path cuppa_resources

    output:
    tuple val(meta), path('*csv'), emit: csv
    path 'versions.yml'          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def has_rna = isofox_dir
    def has_dna = purple_dir || linx_dir
    if (has_dna && (!purple_dir || !linx_dir)) {
        exit 1, "ERROR: CUPPA_CLASSIFIER: DNA classification requires both PURPLE and LINX inputs"
    }

    def categories_val
    if (has_dna && has_rna) {
        categories_val = 'ALL'
    } else if (has_dna) {
        categories_val = 'DNA'
    } else if (has_rna) {
        categories_val = 'RNA'
    } else {
        exit 1, "ERROR: CUPPA_CLASSIFIER: either PURPLE and LINX inputs or Isofox inputs are required"
    }

    """
    # Symlink input files into a single directory
    mkdir -p sample_data/
    if [[ ${categories_val} == 'DNA' || ${categories_val} == 'ALL' ]]; then
        find -L ${purple_dir} ${linx_dir} ${virusinterpreter} -mindepth 1 -maxdepth 1 -type f -exec ln -fs ../{} sample_data/ \\;
    fi

    if [[ ${categories_val} == 'RNA' ]]; then
        find -L ${isofox_dir} -mindepth 1 -maxdepth 1 -type f -exec ln -fs ../{} sample_data/ \\;
    elif [[ ${categories_val} == 'ALL' ]]; then
        # NOTE(SW): CUPPA requires that the WTS sample name matches the WGS sample name
        for fp in \$(find -L ${isofox_dir} -mindepth 1 -maxdepth 1 -type f); do
            fn_out=\$(sed 's/^${meta.id_wts}/${meta.id}/' <<< \${fp##*/});
            cp \${fp} sample_data/\${fn_out}
        done;
        # Rename identifier in the summary file
        sed -i 's/^${meta.id_wts}/${meta.id}/g' sample_data/${meta.id}.isf.summary.csv
    fi;

    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            -categories ${categories_val} \\
            -ref_data_dir ${cuppa_resources} \\
            -sample_data ${meta.id} \\
            -sample_data_dir sample_data/ \\
            -ref_genome_version ${ref_genome_ver} \\
            -output_dir ./

    # NOTE(SW): hard coded since there is no reliable way to obtain version information.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cuppa: 1.8
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.cup.data.csv
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
