process ISOFOX_PANEL_NORMALISATION {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-isofox:1.7.1--hdfd78af_0' :
        'biocontainers/hmftools-isofox:1.7.1--hdfd78af_0' }"

    input:
    path 'isofox_dirs.*'
    val genome_ver
    path gene_ids
    path gene_distribution

    output:
    path 'isofox.gene_normalisation.*.csv', emit: isofox_normalisation
    path 'versions.yml'                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p inputs/
    for fp in \$(find -L isofox_dirs.* -name '*.gene_data.csv'); do ln -sf ../\${fp} inputs/; done

    (
       echo SampleId;
       find inputs/ -name '*csv' | sed 's#^.*/\\(.*\\).isf.gene_data.csv#\\1#';
    ) > sample_ids.txt

    java -cp /usr/local/share/hmftools-isofox-1.7.1-0/isofox.jar \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        com.hartwig.hmftools.isofox.cohort.CohortAnalyser \\
        ${args} \\
        -sample_data_file sample_ids.txt \\
        -root_data_dir inputs/ \\
        -analyses PANEL_TPM_NORMALISATION \\
        -gene_id_file ${gene_ids} \\
        -gene_distribution_file ${gene_distribution} \\
        -output_dir ./ \\
        -log_level ${params.module_log_level}

    mv isofox.panel_gene_normalisation.csv isofox.gene_normalisation.${genome_ver}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isofox: \$(isofox -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch isofox.gene_normalisation.${genome_ver}.csv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
