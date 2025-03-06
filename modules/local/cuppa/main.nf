process CUPPA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-cuppa:2.3.1--py311r42hdfd78af_0' :
        'biocontainers/hmftools-cuppa:2.3.1--py311r42hdfd78af_0' }"

    input:
    tuple val(meta), path(isofox_dir), path(purple_dir), path(linx_dir), path(virusinterpreter_dir)
    val genome_ver
    path cuppa_alt_sj
    path cuppa_classifier
    val categories

    output:
    tuple val(meta), path('cuppa/'), emit: cuppa_dir
    path 'versions.yml'            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def isofox_dir_arg = isofox_dir ? "-isofox_dir isofox_dir__prepared/" : ""
    def ref_alt_sj_sites_arg = isofox_dir ? "-ref_alt_sj_sites ${cuppa_alt_sj}" : ""

    """
    if [[ -n "${isofox_dir}" ]]; then
        # NOTE(SW): CUPPA requires that the RNA sample name matches the DNA sample name
        mkdir -p isofox_dir__prepared/
        for fp in ${isofox_dir}/*; do
            cp -L \${fp} isofox_dir__prepared/\$(sed 's/${meta.sample_rna_id}/${meta.sample_id}/' <<< \${fp##*/});
        done;
    fi;

    mkdir -p cuppa/

    # Extract input features
    cuppa \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        com.hartwig.hmftools.cup.prep.CuppaDataPrep \\
        ${args} \\
        -sample ${meta.sample_id} \\
        -categories ${categories} \\
        -ref_genome_version ${genome_ver} \\
        -output_dir cuppa/ \\
        -purple_dir ${purple_dir} \\
        -linx_dir ${linx_dir} \\
        -virus_dir ${virusinterpreter_dir} \\
        ${isofox_dir_arg} \\
        ${ref_alt_sj_sites_arg}

    # Make predictions
    python -m cuppa.predict \\
        --sample_id ${meta.sample_id} \\
        --classifier_path ${cuppa_classifier} \\
        --features_path cuppa/${meta.sample_id}.cuppa_data.tsv.gz \\
        --output_dir cuppa/ \\
        --clf_group ${categories}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cuppa: \$(cuppa -version | sed -n '/Cuppa version/ { s/^.* //p }')
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
