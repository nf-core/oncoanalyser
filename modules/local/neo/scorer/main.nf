process NEO_SCORER {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/neo:1.1_beta--0'

    input:
    tuple val(meta), path(purple_dir), path(isofox_dir), path(lilac_dir), path(annotate_fusions_dir), path(neo_finder_dir)
    path genome_fasta
    val genome_ver
    path ensembl_data_resources
    path neo_resources, stageAs: 'neo_reference_data'
    path cohort_tpm_medians

    output:
    tuple val(meta), path('neo/'), emit: neo_scorer_dir
    path 'versions.yml'          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def rna_sample_arg = meta.containsKey('sample_rna_id') ? "-rna_sample ${meta.sample_rna_id}" : ''
    def rna_somatic_vcf_arg = meta.containsKey('sample_rna_id') ? "-rna_somatic_vcf ${purple_dir}/${meta.sample_id}.sage_append.vcf.gz" : ''

    // NeoScorer expects the fusion-neoepitopes which Isofox has annotated with RNA to be in the Isofox directory, so put them
    // and the standard Isofox files (just TPM is used) into a new shared directory
    // ie isofox_neo_dir + neo_finder_dir -> new directory for isofox data -> passed into -isofox_dir
    def isofox_dir_arg = meta.containsKey('sample_rna_id') ? "-isofox_dir /path/isofox_combined_dir" : ''

    """
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -cp ${task.ext.jarPath} \\
        com.hartwig.hmftools.neo.score.NeoScorer \\
            ${args} \\
            -sample ${meta.sample_id} \\
            ${rna_sample_arg} \\
            -purple_dir ${purple_dir} \\
            -lilac_dir ${lilac_dir} \\
            ${isofox_dir_arg} \\
            ${rna_somatic_vcf_arg} ]]
            -neo_dir ${neo_finder_dir} \\
            -ref_genome ${genome_fasta} \\
            -ref_genome_version ${genome_ver} \\
            -ensembl_data_dir ${ensembl_data_resources} \\
            -score_file_dir ${neo_resources} \\
            -cancer_tpm_medians_file ${cohort_tpm_medians} \\
            -output_dir ${output_dir} \\
            -log_debug \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        neo: \$(java -jar ${task.ext.jarPath} -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p neo/
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}

