process COBALT_PANEL_NORMALISATION {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-cobalt:3.0--hdfd78af_0' :
        'biocontainers/hmftools-cobalt:3.0--hdfd78af_0' }"

    input:
    tuple path('amber_dir.*'), path('cobalt_dir.*')
    val genome_ver
    path gc_profile
    path copy_number_percentiles
    path target_regions_bed

    output:
    path 'cobalt.region_normalisation.*.tsv'                             , topic: cobalt_normalisation_tsv
    tuple val([:]), val('cobalt_panel_normalisation'), path('.command.*'), topic: command_files
    path 'versions.yml'                                                  , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def wgs_copy_number_percentiles_arg = copy_number_percentiles ? "-wgs_copy_number_percentiles ${copy_number_percentiles}" : ''

    """
    mkdir -p inputs/

    for fp in \$(find -L amber_dir.* cobalt_dir.* -type f ! -name '*.version'); do
        ln -sf ../\${fp} inputs/\${fp##*/};
    done

    (
        echo SampleId
        basename -s .amber.baf.tsv.gz -a inputs/*.amber.baf.tsv.gz
    ) > sample_ids.txt

    cobalt \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        com.hartwig.hmftools.cobalt.norm.NormalisationFileBuilder \\
        ${args} \\
        -sample_id_file sample_ids.txt \\
        -amber_dir inputs/ \\
        -cobalt_dir inputs/ \\
        -ref_genome_version ${genome_ver} \\
        -gc_profile ${gc_profile} \\
        ${wgs_copy_number_percentiles_arg} \\
        -target_regions_bed ${target_regions_bed} \\
        ${log_level_arg} \\
        -output_file cobalt.region_normalisation.${genome_ver}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobalt: \$(cobalt -version | sed -n '/^Cobalt version / { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
        r: \$(R --version | sed -n '/^R version/ { s/^.*version //; s/ .*//p }')
        r-dplyr: \$(Rscript -e 'packageVersion("dplyr") |> as.character() |> writeLines()')
        bioconductor-copynumber: \$(Rscript -e 'packageVersion("copynumber") |> as.character() |> writeLines()')
    END_VERSIONS
    """

    stub:
    """
    touch cobalt.region_normalisation.${genome_ver}.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
