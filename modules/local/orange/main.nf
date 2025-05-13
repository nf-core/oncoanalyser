process ORANGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-orange:3.8.1--hdfd78af_0' :
        'biocontainers/hmftools-orange:3.8.1--hdfd78af_0' }"

    input:
    tuple val(meta),
        path(bamtools_somatic_dir),
        path(bamtools_germline_dir),
        path(sage_somatic_dir),
        path(sage_germline_dir),
        path(smlv_somatic_vcf),
        path(smlv_germline_vcf),
        path(purple_dir),
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
    path cohort_mapping
    path cohort_percentiles
    path known_fusion_data
    path driver_gene_panel
    path ensembl_data_resources
    path sigs_etiology
    path isofox_alt_sj
    path isofox_gene_distribution
    val pipeline_version

    output:
    tuple val(meta), path('output/*.orange.pdf') , emit: pdf, optional: true
    tuple val(meta), path('output/*.orange.json'), emit: json, optional: true
    path 'versions.yml'                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def xmx_mod = task.ext.xmx_mod ?: 0.95

    def pipeline_version_str = pipeline_version ?: 'not specified'

    def run_mode = Utils.getEnumFromString(params.mode, Constants.RunMode);
    def experiment_type = (run_mode === Constants.RunMode.WGTS) ? 'WGS' : 'PANEL'

    def virus_dir_arg = virusinterpreter_dir ? "-virus_dir ${virusinterpreter_dir}" : ''
    def lilac_dir_arg = lilac_dir ? "-lilac_dir ${lilac_dir}" : ''
    def chord_dir_arg = chord_dir ? "-chord_dir ${chord_dir}" : ''
    def sigs_dir_arg = sigs_dir ? "-sigs_dir ${sigs_dir}" : ''
    def cuppa_dir_arg = cuppa_dir ? "-cuppa_dir ${cuppa_dir}" : ''
    def peach_dir_arg = peach_dir ? "-peach_dir ${peach_dir}" : ''
    def plot_dir = linx_somatic_plot_dir.resolve('reportable/').toUriString().replaceAll('/$', '')

    def tumor_metrics_arg = "-tumor_metrics_dir ${bamtools_somatic_dir}"
    def normal_metrics_arg = bamtools_germline_dir ? "-ref_metrics_dir ${bamtools_germline_dir}" : ''

    def normal_id_arg = meta.containsKey('normal_dna_id') ? "-reference_sample_id ${meta.normal_dna_id}" : ''
    def normal_sage_dir = sage_germline_dir ? "-sage_germline_dir ${sage_germline_dir}" : ''
    def normal_linx_arg = linx_germline_anno_dir ? "-linx_germline_dir ${linx_germline_anno_dir}" : ''

    def rna_id_arg = meta.containsKey('tumor_rna_id') ? "-rna_sample_id ${meta.tumor_rna_id}" : ''
    def isofox_dir_arg = isofox_dir ? '-isofox_dir isofox_dir__prepared/' : ''

    def isofox_gene_distribution_arg = isofox_gene_distribution ? "-isofox_gene_distribution ${isofox_gene_distribution}" : ''
    def isofox_alt_sj_arg = isofox_alt_sj ? "-isofox_alt_sj_cohort ${isofox_alt_sj}" : ''

    // NOTE(SW): DOID label: 162 [cancer]; Hartwig cohort group: unknown
    def doid_arg = meta.cancer_type ?: '162'

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
    if [[ ! -e ${plot_dir}/ ]]; then
        mkdir -p ${plot_dir}/;
    fi;

    mkdir -p output/

    orange \\
        -Xmx${Math.round(task.memory.bytes * xmx_mod)} \\
        ${args} \\
        \\
        -add_disclaimer \\
        -pipeline_version_file pipeline_version.txt \\
        -experiment_type ${experiment_type} \\
        \\
        -tumor_sample_id ${meta.tumor_id} \\
        -primary_tumor_doids ${doid_arg} \\
        -sage_dir ${sage_somatic_dir} \\
        -purple_dir \${purple_dir_local} \\
        -purple_plot_dir \${purple_dir_local}/plot/ \\
        -linx_dir ${linx_somatic_anno_dir} \\
        -linx_plot_dir ${plot_dir}/ \\
        ${virus_dir_arg} \\
        ${lilac_dir_arg} \\
        ${chord_dir_arg} \\
        ${sigs_dir_arg} \\
        ${cuppa_dir_arg} \\
        ${peach_dir_arg} \\
        \\
        ${normal_id_arg} \\
        ${normal_metrics_arg} \\
        ${tumor_metrics_arg} \\
        ${normal_sage_dir} \\
        ${normal_linx_arg} \\
        \\
        ${rna_id_arg} \\
        ${isofox_dir_arg} \\
        \\
        -ref_genome_version ${genome_ver} \\
        -doid_json ${disease_ontology} \\
        -cohort_mapping_tsv ${cohort_mapping} \\
        -cohort_percentiles_tsv ${cohort_percentiles} \\
        -known_fusion_file ${known_fusion_data} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -signatures_etiology_tsv ${sigs_etiology} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        ${isofox_gene_distribution_arg} \\
        ${isofox_alt_sj_arg} \\
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
