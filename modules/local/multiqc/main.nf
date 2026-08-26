process MULTIQC {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d8bcaa0f18f0068460deb4841052ef5429108a27:f0d5b44eec0233d9bc8357f24755076bceef041a-0' :
        'biocontainers/mulled-v2-d8bcaa0f18f0068460deb4841052ef5429108a27:f0d5b44eec0233d9bc8357f24755076bceef041a-0' }"

    input:
    tuple val(meta), val(fs_group_ids), path(fs, name: 'sample_flat/*/*')
    path multiqc_files_other, stageAs: "other/*"
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)

    output:
    path '*multiqc_report.html'                        , topic: multiqc_report
    path '*_data'                                      , topic: multiqc_data
    path '*_plots'                                     , topic: multiqc_plots, optional:true
    tuple val(meta), val('multiqc'), path('.command.*'), topic: command_files
    // NOTE(SW): cannot published otherwise produce ~ consume on topic channel becomes circular and blocks
    //path "versions.yml"                                , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def filename_arg = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : ''

    def config_arg = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config_arg = extra_multiqc_config ? "--config $extra_multiqc_config" : ''

    def logo_arg = multiqc_logo ? "--cl-config 'custom_logo: \"${multiqc_logo}\"'" : ''

    // Prepare sample linkage for process
    def group_id_data = meta.collectEntries { gid, _sns -> [gid, "sample/${gid}"] }

    def group_id_mapping = [fs_group_ids, fs]
        .transpose()
        .collect { gid, gfn ->
            def fd = group_id_data[gid]
            "${gfn},${fd}"
        }
        .join('\n')

    def manifest = meta
        .collect { gid, sns ->
            def fd = group_id_data[gid]
            "${gid},${sns.normal_dna_id},${sns.tumor_dna_id},${sns.tumor_rna_id},${fd}"
        }
        .join('\n')

    """
    # Write manifest to disk and structure sample inputs
    echo "${manifest}" > manifest.csv

    while IFS=',' read -r fp fd; do
        mkdir -p \${fd};
        ln -s ../../\${fp%/} \${fd};
    done <<< "${group_id_mapping}"

    # Run MultiQC
    multiqc \\
        ${args} \\
        --force \\
        ${config_arg} \\
        ${extra_config_arg} \\
        ${logo_arg} \\
        ${filename_arg} \\
        --module hmftools_amber \\
        --module hmftools_bamtools \\
        --module hmftools_purple \\
        --module hmftools_rna \\
        --case_information_fp manifest.csv \\
        other/ \\
        sample/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p multiqc_data/
    mkdir -p multiqc_plots/

    touch multiqc_report.html
    touch multiqc_data/.stub
    touch multiqc_plots/.stub

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
