process COBALT_PANEL_NORMALISATION {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-cobalt:2.0--hdfd78af_0' :
        'biocontainers/hmftools-cobalt:2.0--hdfd78af_0' }"

    input:
    tuple path('amber_dir.*'), path('cobalt_dir.*')
    val genome_ver
    path gc_profile
    path target_regions_bed

    output:
    path 'cobalt.region_normalisation.*.tsv', emit: cobalt_normalisation
    path 'versions.yml'                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p inputs/

    for fp in \$(find -L amber_dir.* cobalt_dir.* -type f ! -name '*.version'); do
        ln -s ../\${fp} inputs/\${fp##*/};
    done

    (
        echo SampleId;
        find -L inputs/ -type f -name '*.amber.baf.tsv.gz' | sed 's#inputs/##; s#\\.amber\\..*\$##; s#\\.cobalt\\..*\$##' | sort -V | uniq;
    ) > sample_ids.txt

    java -cp /usr/local/share/hmftools-cobalt-2.0-0/cobalt.jar \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        com.hartwig.hmftools.cobalt.norm.NormalisationFileBuilder \\
            ${args} \\
            -sample_id_file sample_ids.txt \\
            -amber_dir inputs/ \\
            -cobalt_dir inputs/ \\
            -ref_genome_version ${genome_ver} \\
            -gc_profile ${gc_profile} \\
            -target_regions_bed ${target_regions_bed} \\
            -output_file cobalt.region_normalisation.${genome_ver}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobalt_panel_normalisation: \$(cobalt -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    touch cobalt.region_normalisation.${genome_ver}.tsv

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
