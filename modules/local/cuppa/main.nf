process CUPPA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-cuppa:2.5.1--py311r42hdfd78af_0' :
        'biocontainers/hmftools-cuppa:2.5.1--py311r42hdfd78af_0' }"

    input:
    tuple val(meta), path(isofox_dir), path(purple_dir), path(linx_dir), path(virusinterpreter_dir)
    val genome_ver
    path cuppa_alt_sj
    path cuppa_classifier
    val categories

    output:
    tuple val(meta), path('cuppa/'), emit: cuppa_dir
    path 'versions.yml'            , emit: versions
    path '.command.*'              , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.75

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def purple_dir_arg = purple_dir ? "-purple_dir ${purple_dir}" : ''
    def linx_dir_arg = linx_dir ? "-linx_dir ${linx_dir}" : ''
    def virusinterpreter_dir_arg = virusinterpreter_dir ? "-virus_dir ${virusinterpreter_dir}" : ''

    def isofox_dir_arg = ''
    def rna_sample_arg = ''
    def ref_alt_sj_sites_arg = ''

    if (isofox_dir) {
        isofox_dir_arg = "-isofox_dir ${isofox_dir}"
        rna_sample_arg = "-rna_sample ${meta.sample_rna_id}"
        ref_alt_sj_sites_arg = "-ref_alt_sj_sites ${cuppa_alt_sj}"
    }

    """
    mkdir -p cuppa/

    # Extract input features
    cuppa \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        com.hartwig.hmftools.cup.prep.CuppaDataPrep \\
        ${args} \\
        -sample ${meta.sample_id} \\
        ${rna_sample_arg} \\
        -categories ${categories} \\
        ${purple_dir_arg} \\
        ${linx_dir_arg} \\
        ${virusinterpreter_dir_arg} \\
        ${isofox_dir_arg} \\
        -ref_genome_version ${genome_ver} \\
        ${ref_alt_sj_sites_arg} \\
        ${log_level_arg} \\
        -output_dir cuppa/

    # Make predictions
    python -m cuppa.predict \\
        ${args2} \\
        --sample_id ${meta.sample_id} \\
        --features_path cuppa/${meta.sample_id}.cuppa_data.tsv.gz \\
        --clf_group ${categories} \\
        --classifier_path ${cuppa_classifier} \\
        --output_dir cuppa/

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
