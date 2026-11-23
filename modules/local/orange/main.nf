process ORANGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-orange:5.0.3--hdfd78af_0' :
        'biocontainers/hmftools-orange:5.0.3--hdfd78af_0' }"

    input:
    tuple val(meta),
        path(sage_dir_somatic, stageAs: 'sage_somatic'),
        path(sage_dir_germline, stageAs: 'sage_germline'),
        path(sage_plot_dir_somatic),
        path(sage_append_dir),
        path(purple_dir),
        path(qsee_dir),
        path(linx_annotation_dir_somatic),
        path(linx_plot_dir_reportable_somatic),
        path(linx_annotation_dir_germline),
        path(viridian_dir),
        path(chord_dir),
        path(sigs_dir),
        path(lilac_dir),
        path(cuppa_dir),
        path(peach_dir),
        path(isofox_dir)
    val genome_ver
    path disease_ontology
    val pipeline_version
    val sequencing_platform
    val targeted_mode
    val panel

    output:
    tuple val(meta), path('output/*.orange.pdf')      , topic: orange_pdf, optional: true
    tuple val(meta), path('output/*.orange.json')     , topic: orange_json, optional: true
    tuple val(meta), val('orange'), path('.command.*'), topic: command_files
    path 'versions.yml'                               , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''

    def pipeline_version_str = pipeline_version ?: 'not specified'

    def experiment_type = 'WGS'
    def panel_name_arg = ''
    if (targeted_mode) {
        experiment_type = 'PANEL'
        panel_name_arg = "-panel_name ${panel.toUpperCase()}"
    }

    def primary_tumor_location_arg = meta.cancer_type ? "-primary_tumor_location ${meta.cancer_type}" : ''

    def reference_arg = meta.containsKey('normal_dna_id') ? "-reference ${meta.normal_dna_id}" : ''
    def sage_germline_dir_arg = sage_dir_germline ? "-sage_germline_dir ${sage_dir_germline}" : ''
    def linx_germline_dir_arg = linx_annotation_dir_germline ? "-linx_germline_dir ${linx_annotation_dir_germline}" : ''

    def sage_plot_dir_arg = sage_plot_dir_somatic ? "-sage_plot_dir ${sage_plot_dir_somatic}" : ''
    def viridian_dir_arg = viridian_dir ? "-virus_dir ${viridian_dir}" : ''
    def lilac_dir_arg = lilac_dir ? "-lilac_dir ${lilac_dir}" : ''
    def chord_dir_arg = chord_dir ? "-chord_dir ${chord_dir}" : ''
    def sigs_dir_arg = sigs_dir ? "-sigs_dir ${sigs_dir}" : ''
    def cuppa_dir_arg = cuppa_dir ? "-cuppa_dir ${cuppa_dir}" : ''
    def peach_dir_arg = peach_dir ? "-peach_dir ${peach_dir}" : ''

    def rna_sample_arg = meta.containsKey('tumor_rna_id') ? "-rna_sample ${meta.tumor_rna_id}" : ''
    def isofox_dir_arg = isofox_dir ? "-isofox_dir ${isofox_dir}" : ''
    def purple_rna_append_vcf_arg = sage_append_dir ? "-purple_rna_append_vcf ${sage_append_dir}/${meta.tumor_id}.sage.append.vcf.gz" : ''

    """
    echo "${pipeline_version_str}" > pipeline_version.txt

    # Set input plot directory and create it doesn't exist. See the LINX visualiser module for further info.
    if [[ ! -e ${linx_plot_dir_reportable_somatic}/ ]]; then
        mkdir -p ${linx_plot_dir_reportable_somatic}/;
    fi;

    mkdir -p output/

    orange \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        \\
        -add_disclaimer \\
        -pipeline_version_file pipeline_version.txt \\
        -experiment_type ${experiment_type} \\
        -sequencing_type ${sequencing_platform.toUpperCase()} \\
        ${panel_name_arg} \\
        ${primary_tumor_location_arg} \\
        \\
        -tumor ${meta.tumor_id} \\
        -sage_dir ${sage_dir_somatic} \\
        -purple_dir ${purple_dir} \\
        -purple_plot_dir ${purple_dir}/plot/ \\
        -qsee_dir ${qsee_dir} \\
        -linx_dir ${linx_annotation_dir_somatic} \\
        -linx_plot_dir ${linx_plot_dir_reportable_somatic}/ \\
        \\
        ${reference_arg} \\
        ${sage_germline_dir_arg} \\
        ${linx_germline_dir_arg} \\
        \\
        ${sage_plot_dir_arg} \\
        ${viridian_dir_arg} \\
        ${lilac_dir_arg} \\
        ${chord_dir_arg} \\
        ${sigs_dir_arg} \\
        ${cuppa_dir_arg} \\
        ${peach_dir_arg} \\
        \\
        ${rna_sample_arg} \\
        ${isofox_dir_arg} \\
        ${purple_rna_append_vcf_arg} \\
        \\
        -ref_genome_version ${genome_ver} \\
        -doid_json ${disease_ontology} \\
        ${log_level_arg} \\
        -output_dir output/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orange: \$(orange -version | sed -n '/^Orange version / { s/^.* //p }')
        java: \$(java --version | sed -n '/^openjdk/ { s/^.*openjdk //; s/ .*//p }')
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
