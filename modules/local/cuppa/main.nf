process CUPPA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-cuppa:1.8.1--hdfd78af_0' :
        'biocontainers/hmftools-cuppa:1.8.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(isofox_dir), path(purple_dir), path(linx_dir), path(virusinterpreter_dir)
    val genome_ver
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

    elif [[ ${classifier} == 'ALL' ]]; then
        # NOTE(SW): CUPPA requires that the RNA sample name matches the DNA sample name
        for fp in \$(find -L ${isofox_dir} -maxdepth 1 -type f); do
            fn_out=\$(sed 's/^${meta.sample_rna_id}/${meta.sample_id}/' <<< \${fp##*/});
            cp \${fp} sample_data/\${fn_out}
        done;
        # Rename identifier in the summary file
        sed -i 's/^${meta.sample_rna_id}/${meta.sample_id}/g' sample_data/${meta.sample_id}.isf.summary.csv;
    fi;

    # Use symlink to remove genome version suffix (e.g. cuppa_classifier.37.pickle.gz -> cuppa_classifier.pickle.gz)
    ln -sf \$(find -L ${cuppa_resources} -type f -name 'cuppa_classifier*.pickle.gz') ${cuppa_resources}/cuppa_classifier.pickle.gz
    ln -sf \$(find -L ${cuppa_resources} -type f -name 'alt_sj.selected_loci*.tsv.gz') ${cuppa_resources}/alt_sj.selected_loci.tsv.gz

    mkdir -p cuppa/

    # Extract input features
    cuppa \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -categories ${classifier} \\
        -ref_genome_version ${genome_ver} \\
        -sample_data_dir sample_data/ \\
        -output_dir cuppa/ \\
        -ref_alt_sj_sites "${cuppa_resources}/alt_sj.selected_loci.tsv.gz"

    # Make predictions
    pyenv activate pycuppa_env
    python3 -m cuppa.predict \\
        --sample_id ${meta.sample_id} \\
        --classifier_path ${cuppa_resources}/cuppa_classifier.pickle.gz \\
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
