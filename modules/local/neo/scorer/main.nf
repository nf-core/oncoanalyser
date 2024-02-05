process NEO_SCORER {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/neo:1.2_beta--1'

    input:
    tuple val(meta), path(isofox_dir), path(purple_dir), path(sage_vcf), path(lilac_dir), path(neo_finder_dir), path(annotate_fusions)
    path ensembl_data_resources
    path neo_resources, stageAs: 'neo_reference_data'
    path cohort_tpm_medians

    output:
    tuple val(meta), path('neo_scorer/'), emit: neo_scorer_dir
    path 'versions.yml'                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def rna_sample_arg = meta.containsKey('sample_rna_id') ? "-rna_sample ${meta.sample_rna_id}" : ''
    def rna_somatic_vcf_arg = meta.containsKey('sample_rna_id') ? "-rna_somatic_vcf ${sage_vcf}" : ''

    """
    isofox_dir_arg=''
    if [[ -n "${isofox_dir}" ]]; then
        isofox_dir_local=isofox__prepared/;

        cp -rL ${isofox_dir} \${isofox_dir_local}/;
        cp -r ${annotate_fusions} \${isofox_dir_local}/;

        isofox_dir_arg="-isofox_dir \${isofox_dir_local}";
    fi;

    mkdir -p neo_scorer/

    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -cp ${task.ext.jarPath} \\
        com.hartwig.hmftools.neo.score.NeoScorer \\
            ${args} \\
            -sample ${meta.sample_id} \\
            ${rna_sample_arg} \\
            \${isofox_dir_arg} \\
            -purple_dir ${purple_dir} \\
            ${rna_somatic_vcf_arg} \\
            -lilac_dir ${lilac_dir} \\
            -neo_dir ${neo_finder_dir} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -score_file_dir ${neo_resources} \\
            -cancer_tpm_medians_file ${cohort_tpm_medians} \\
            -log_debug \\
            -output_dir neo_scorer/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        neo: \$(java -jar ${task.ext.jarPath} -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p neo_scorer/
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}

