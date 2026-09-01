process ISOFOX_PANEL_NORMALISATION {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-isofox:2.0.1--hdfd78af_0' :
        'biocontainers/hmftools-isofox:2.0.1--hdfd78af_0' }"

    input:
    path 'isofox_dirs.*'
    val genome_ver
    path gene_ids
    path gene_distribution

    output:
    path 'isofox.gene_normalisation.*.csv'                               , topic: isofox_normalisation_csv
    tuple val([:]), val('isofox_panel_normalisation'), path('.command.*'), topic: command_files
    path 'versions.yml'                                                  , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    """
    mkdir -p inputs/
    for fp in \$(find -L isofox_dirs.* -name '*.isf.gene_data.tsv'); do ln -sf ../\${fp} inputs/; done

    (
       echo SampleId
       basename -s .isf.gene_data.tsv -a inputs/*.isf.gene_data.tsv
    ) > sample_ids.txt

    isofox \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        com.hartwig.hmftools.isofox.cohort.CohortAnalyser \\
        ${args} \\
        -sample_data_file sample_ids.txt \\
        -root_data_dir inputs/ \\
        -analyses PANEL_TPM_NORMALISATION \\
        -gene_id_file ${gene_ids} \\
        -gene_distribution_file ${gene_distribution} \\
        ${log_level_arg} \\
        -output_dir ./

    mv isofox.panel_gene_normalisation.csv isofox.gene_normalisation.${genome_ver}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isofox: \$(isofox -version | sed -n '/^Isofox version / { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
    END_VERSIONS
    """

    stub:
    """
    touch isofox.gene_normalisation.${genome_ver}.csv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
