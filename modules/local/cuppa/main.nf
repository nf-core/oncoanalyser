process CUPPA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-cuppa:1.8.1--hdfd78af_0' :
        'quay.io/biocontainers/hmftools-cuppa:1.8.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(isofox_dir), path(purple_dir), path(linx_dir), path(virusinterpreter_dir)
    val ref_genome_ver
    path cuppa_resources, stageAs: 'cuppa_reference_data'
    val classifier

    output:
    tuple val(meta), path('cuppa/'), emit: cuppa_dir
    path 'versions.yml'            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    # Symlink input files into a single directory
    mkdir -p sample_data/
    if [[ ${classifier} == 'DNA' || ${classifier} == 'ALL' ]]; then
        find -L ${purple_dir} ${linx_dir} ${virusinterpreter_dir} -maxdepth 1 -type f -exec ln -fs ../{} sample_data/ \\;
    fi

    if [[ ${classifier} == 'RNA' ]]; then
        find -L ${isofox_dir} -maxdepth 1 -type f -exec ln -fs ../{} sample_data/ \\;
    elif [[ ${classifier} == 'ALL' ]]; then
        # NOTE(SW): CUPPA requires that the RNA sample name matches the DNA sample name
        for fp in \$(find -L ${isofox_dir} -maxdepth 1 -type f); do
            fn_out=\$(sed 's/^${meta.sample_rna_id}/${meta.sample_id}/' <<< \${fp##*/});
            cp \${fp} sample_data/\${fn_out}
        done;
        # Rename identifier in the summary file
        sed -i 's/^${meta.sample_rna_id}/${meta.sample_id}/g' sample_data/${meta.sample_id}.isf.summary.csv;
    fi;

    mkdir -p cuppa/

    cuppa \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -sample ${meta.sample_id} \\
        -sample_data_dir sample_data/ \\
        -categories ${classifier} \\
        -ref_data_dir ${cuppa_resources} \\
        -ref_genome_version ${ref_genome_ver} \\
        -create_pdf \\
        -output_dir cuppa/

    if [[ ${classifier} == 'DNA' || ${classifier} == 'ALL' ]]; then
        cuppa-chart \\
            -sample ${meta.sample_id} \\
            -sample_data cuppa/${meta.sample_id}.cup.data.csv \\
            -output_dir cuppa/;
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cuppa: \$(cuppa | sed -n '1s/^.* //p')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p cuppa/
    touch cuppa/${meta.sample_id}.cup.data.csv
    touch cuppa/${meta.sample_id}.cuppa.conclusion.txt
    touch cuppa/${meta.sample_id}_cup_report.pdf
    touch cuppa/${meta.sample_id}.cup.report.summary.png
    touch cuppa/${meta.sample_id}.cup.report.features.png
    touch cuppa/${meta.sample_id}.cuppa.chart.png
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
