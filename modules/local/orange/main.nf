process ORANGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-orange:4.1.2--hdfd78af_0' :
        'biocontainers/hmftools-orange:4.1.2--hdfd78af_0' }"

    input:
    tuple val(meta),
        path(sage_somatic_dir, stageAs: 'sage_somatic'),
        path(sage_germline_dir, stageAs: 'sage_germline'),
        path(smlv_somatic_vcf),
        path(smlv_germline_vcf),
        path(purple_dir),
        path(qsee_dir),
        path(linx_somatic_anno_dir),
        path(linx_somatic_plot_dir),
        path(linx_germline_anno_dir),
        path(virusinterpreter_dir),
        path(chord_dir),
        path(sigs_dir),
        path(lilac_dir),
        path(cuppa_dir),
        path(peach_dir),
        path(isofox_dir)
    val genome_ver
    path disease_ontology
    val pipeline_version
    val sequencing_type
    val targeted_mode

    output:
    tuple val(meta), path('output/*.orange.pdf') , emit: pdf, optional: true
    tuple val(meta), path('output/*.orange.json'), emit: json, optional: true
    path 'versions.yml'                          , emit: versions
    path '.command.*'                            , emit: command_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def pipeline_version_str = pipeline_version ?: 'not specified'
    def experiment_type = targeted_mode ? 'PANEL' : 'WGS'
    def cancer_type_arg = meta.cancer_type ? "-primary_tumor_location ${meta.cancer_type}" : ''

    // Tumor sample
    def linx_plot_dir = linx_somatic_plot_dir.resolve('reportable/').toUriString().replaceAll('/$', '')

    // Normal sample
    def normal_id_arg = meta.containsKey('normal_dna_id') ? "-reference ${meta.normal_dna_id}" : ''
    def normal_sage_dir_arg = sage_germline_dir ? "-sage_germline_dir ${sage_germline_dir}" : ''
    def normal_linx_arg = linx_germline_anno_dir ? "-linx_germline_dir ${linx_germline_anno_dir}" : ''

    // Optional tools
    def virus_dir_arg = virusinterpreter_dir ? "-virus_dir ${virusinterpreter_dir}" : ''
    def lilac_dir_arg = lilac_dir ? "-lilac_dir ${lilac_dir}" : ''
    def chord_dir_arg = chord_dir ? "-chord_dir ${chord_dir}" : ''
    def sigs_dir_arg = sigs_dir ? "-sigs_dir ${sigs_dir}" : ''
    def cuppa_dir_arg = cuppa_dir ? "-cuppa_dir ${cuppa_dir}" : ''
    def peach_dir_arg = peach_dir ? "-peach_dir ${peach_dir}" : ''

    // RNA sample
    def rna_id_arg = meta.containsKey('tumor_rna_id') ? "-rna_sample_id ${meta.tumor_rna_id}" : ''
    def isofox_dir_arg = isofox_dir ? '-isofox_dir isofox_dir__prepared/' : ''

    """
    echo "${pipeline_version_str}" > pipeline_version.txt

    # When WTS data is present, ORANGE expects the somatic SAGE VCF to have appended WTS data; CS indicates this should
    # occur after PURPLE. Since ORANGE only collects the somatic SAGE VCF from the PURPLE output directory, we must
    # prepare accordingly

    # Isofox inputs are also expected to have the tumor sample ID in the filename

    # NOTES(SW): Use of symlinks was causing reliability issues on HPC with Singularity, switched to full file copy instead

    purple_dir_local=${purple_dir}
    if [[ -n "${rna_id_arg}" ]]; then

        purple_dir_local=purple__prepared;

        if [[ -d \${purple_dir_local}/ ]]; then
            rm -r \${purple_dir_local}/;
        fi

        cp -rL ${purple_dir} \${purple_dir_local}/
        cp -L ${smlv_somatic_vcf} \${purple_dir_local}/${meta.tumor_id}.purple.somatic.vcf.gz;

        if [[ -n "${smlv_germline_vcf}" ]]; then
            cp -L ${smlv_germline_vcf} \${purple_dir_local}/${meta.tumor_id}.purple.germline.vcf.gz;
        fi;

        mkdir -p isofox_dir__prepared/;
        for fp in ${isofox_dir}/*; do
            cp -L \${fp} isofox_dir__prepared/\$(sed 's/${meta.tumor_rna_id}/${meta.tumor_id}/' <<< \${fp##*/});
        done;

    fi

    # Set input plot directory and create it doesn't exist. See the LINX visualiser module for further info.
    if [[ ! -e ${linx_plot_dir}/ ]]; then
        mkdir -p ${linx_plot_dir}/;
    fi;

    mkdir -p output/

    orange \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        \\
        -add_disclaimer \\
        -pipeline_version_file pipeline_version.txt \\
        -experiment_type ${experiment_type} \\
        -sequencing_type ${sequencing_type} \\
        ${cancer_type_arg} \\
        \\
        -tumor ${meta.tumor_id} \\
        -sage_dir ${sage_somatic_dir} \\
        -purple_dir \${purple_dir_local} \\
        -purple_plot_dir \${purple_dir_local}/plot/ \\
        -qsee_dir ${qsee_dir} \\
        -linx_dir ${linx_somatic_anno_dir} \\
        -linx_plot_dir ${linx_plot_dir}/ \\
        \\
        ${normal_id_arg} \\
        ${normal_sage_dir_arg} \\
        ${normal_linx_arg} \\
        \\
        ${virus_dir_arg} \\
        ${lilac_dir_arg} \\
        ${chord_dir_arg} \\
        ${sigs_dir_arg} \\
        ${cuppa_dir_arg} \\
        ${peach_dir_arg} \\
        \\
        ${rna_id_arg} \\
        ${isofox_dir_arg} \\
        \\
        -ref_genome_version ${genome_ver} \\
        -doid_json ${disease_ontology} \\
        ${log_level_arg} \\
        -output_dir output/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orange: \$(orange -version | sed -n '/^Orange version / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p output/

    touch output/${meta.tumor_id}.orange.json
    touch output/${meta.tumor_id}.orange.pdf

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
