process NEO_SCORER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-neo:1.3--hdfd78af_0' :
        'biocontainers/hmftools-neo:1.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(isofox_dir), path(purple_dir), path(sage_vcf), path(lilac_dir), path(neo_finder_dir), path(annotated_fusions)
    path ensembl_data_resources
    path neo_resources, stageAs: 'neo_reference_data'
    path cohort_tpm_medians

    output:
    tuple val(meta), path('neo_scorer/')                  , topic: neo_scorer_dir
    tuple val(meta), val('neo_scorer'), path('.command.*'), topic: command_files
    path 'versions.yml'                                   , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def rna_sample_arg = meta.containsKey('sample_rna_id') ? "-rna_sample ${meta.sample_id}" : ''
    def rna_somatic_vcf_arg = meta.containsKey('sample_rna_id') ? "-rna_somatic_vcf ${sage_vcf}" : ''

    def cancer_type_arg = meta.cancer_type ? "-cancer_type ${meta.cancer_type}" : ''

    def isofox_dir_local = 'isofox__prepared/'
    def isofox_dir_arg = isofox_dir ? "-isofox_dir ${isofox_dir_local}" : ''

    """
    mkdir -p neo_scorer/

    # Prepare ISOFOX results
    if [[ -n "${isofox_dir_arg}" ]]; then
      # When provided existing results generated in a RNA-only analysis we must adjust identifier
      if [[ -n "\$(find -L ${isofox_dir} -name '${meta.sample_rna_id}*')" ]]; then
          mkdir -p ${isofox_dir_local}/;
          for e in \$(find -L ${isofox_dir}/*); do
             s=\$(sed 's/^${meta.sample_rna_id}//' <<< \${e##*/});
             ln -s ../\${e} ${isofox_dir_local}/${meta.sample_id}\${s};
          done;
      else
        cp -rL ${isofox_dir} ${isofox_dir_local};
      fi;
      cp -r ${annotated_fusions} ${isofox_dir_local};
    fi

    neo \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        com.hartwig.hmftools.neo.score.NeoScorer \\
        ${args} \\
        -sample ${meta.sample_id} \\
        ${cancer_type_arg} \\
        -purple_dir ${purple_dir} \\
        ${rna_sample_arg} \\
        ${isofox_dir_arg} \\
        ${rna_somatic_vcf_arg} \\
        -lilac_dir ${lilac_dir} \\
        -neo_dir ${neo_finder_dir} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -score_file_dir ${neo_resources} \\
        -cancer_tpm_medians_file ${cohort_tpm_medians} \\
        ${log_level_arg} \\
        -output_dir neo_scorer/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        neo: \$(neo -version | sed -n '/^Neo version / { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p neo_scorer/

    touch neo_scorer/.stub

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
