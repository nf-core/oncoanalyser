process CUPPA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-cuppa:2.3.0_beta--py311r42hdfd78af_0' :
        'biocontainers/hmftools-cuppa:2.3.0_beta--py311r42hdfd78af_0' }"

    input:
    tuple val(meta), path(isofox_dir), path(purple_dir), path(linx_dir), path(virusinterpreter_dir)
    val genome_ver
    path cuppa_alt_sj
    path cuppa_classifier
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

    for file_path in \$(find -L ${purple_dir} ${linx_dir} ${virusinterpreter_dir} -maxdepth 1 -type f -exec realpath {} \\;); do
        ln -sf \${file_path} sample_data/\$(basename \${file_path})
    done;

    if [ ${classifier} == 'ALL' ]; then
        # NOTE(SW): CUPPA requires that the RNA sample name matches the DNA sample name
        for file_path in \$(find -L ${isofox_dir} -maxdepth 1 -type f -exec realpath {} \\;); do
            new_file_name=\$(basename \${file_path} | sed 's/^${meta.sample_rna_id}/${meta.sample_id}/')
            ln -sf \${file_path} sample_data/\${new_file_name}
        done;
    fi;

    mkdir -p cuppa/

    # Extract input features
    cuppa \\
        -Xmx${Math.round(task.memory.bytes * 0.75)} \\
        com.hartwig.hmftools.cup.prep.CuppaDataPrep \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -categories ${classifier} \\
        -ref_genome_version ${genome_ver} \\
        -sample_data_dir sample_data/ \\
        -output_dir cuppa/ \\
        -ref_alt_sj_sites ${cuppa_alt_sj}

    # Make predictions
    python -m cuppa.predict \\
        --sample_id ${meta.sample_id} \\
        --classifier_path ${cuppa_classifier} \\
        --features_path cuppa/${meta.sample_id}.cuppa_data.tsv.gz \\
        --output_dir cuppa/ \\
        --clf_group ${classifier}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cuppa: \$(cuppa | sed -n '1s/^.* //p')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p cuppa/

    touch cuppa/${meta.sample_id}.cuppa_data.tsv.gz
    touch cuppa/${meta.sample_id}.cuppa.pred_summ.tsv
    touch cuppa/${meta.sample_id}.cuppa.vis_data.tsv
    touch cuppa/${meta.sample_id}.cuppa.vis.png

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
