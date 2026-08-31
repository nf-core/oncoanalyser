nextflow.enable.types = true

process MULTIQC {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d8bcaa0f18f0068460deb4841052ef5429108a27:f0d5b44eec0233d9bc8357f24755076bceef041a-0' :
        'biocontainers/mulled-v2-d8bcaa0f18f0068460deb4841052ef5429108a27:f0d5b44eec0233d9bc8357f24755076bceef041a-0' }"

    input:
    tuple(
        case_information: Map,
        fs_case_ids: List<String>,
        fs: List<Path>,
    )
    multiqc_files_other: List<Path>
    multiqc_config: List<Path>
    extra_multiqc_config: List<Path>
    multiqc_logo: List<Path>

    stage:
    stageAs fs, 'sample_flat/*/*'
    stageAs multiqc_files_other, 'other/*'

    topic:
    file('*multiqc_report.html')            >> 'multiqc_report'
    file('*_data')                          >> 'multiqc_data'
    file('*_plots', optional: true)         >> 'multiqc_plots'
    tuple([:], 'multiqc', files('.command.*')) >> 'command_files'
    // NOTE(SW): cannot published otherwise produce ~ consume on topic channel becomes circular and blocks
    // file('versions.yml') >> 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def filename_arg = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : ''

    def config_arg = multiqc_config ? "--config ${multiqc_config.first()}" : ''
    def extra_config_arg = extra_multiqc_config ? "--config ${extra_multiqc_config.first()}" : ''

    def logo_arg = multiqc_logo ? "--cl-config 'custom_logo: \"${multiqc_logo.first()}\"'" : ''

    // Prepare sample linkage for process
    def case_id_data = case_information.collectEntries { gid, _sns -> [gid, "sample/${gid}"] }

    def case_id_mapping = [fs_case_ids, fs]
        .transpose()
        .collect { gid, gfn ->
            def fd = case_id_data[gid]
            "${gfn},${fd}"
        }
        .join('\n')

    def manifest = case_information
        .collect { gid, sns ->
            def fd = case_id_data[gid]
            "${gid},${sns.normal_dna_id},${sns.tumor_dna_id},${sns.tumor_rna_id},${fd}"
        }
        .join('\n')

    """
    # Write manifest to disk and structure sample inputs
    echo "${manifest}" > manifest.csv

    while IFS=',' read -r fp fd; do
        mkdir -p \${fd};
        ln -s ../../\${fp%/} \${fd};
    done <<< "${case_id_mapping}"

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
