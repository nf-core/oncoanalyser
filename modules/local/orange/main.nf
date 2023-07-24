process ORANGE {
    tag "${meta.id}"
    label 'process_single'

    container 'docker.io/scwatts/orange:2.3--0'

    input:
    tuple val(meta), path(bam_metrics_somatic), path(bam_metrics_germline), path(flagstat_somatic), path(flagstat_germline), path(sage_somatic_bqr), path(sage_germline_bqr), path(sage_germline_coverage), path(purple_dir), path(smlv_somatic_vcf), path(smlv_germline_vcf), path(linx_somatic_anno_dir), path(linx_somatic_plot_dir), path(linx_germline_anno_dir), path(virusinterpreter), path(chord_prediction), path(sigs_dir), path(lilac_dir), path(cuppa_dir), path(isofox_dir)
    val genome_ver
    path disease_ontology
    path cohort_mapping
    path cohort_percentiles
    path known_fusion_data
    path driver_gene_panel
    path ensembl_data_resources
    path isofox_alt_sj
    path isofox_gene_distribution
    val pipeline_version

    output:
    tuple val(meta), path('output/*.orange.pdf') , emit: pdf
    tuple val(meta), path('output/*.orange.json'), emit: json
    path 'versions.yml'                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def pipeline_version_str = pipeline_version ?: 'not specified'

    def virusinterpreter_arg = virusinterpreter ? "-annotated_virus_tsv ${virusinterpreter}" : ''
    def chord_arg = chord_prediction ? "-chord_prediction_txt ${chord_prediction}" : ''
    def sigs_arg = sigs_dir ? "-sigs_allocation_tsv ${sigs_dir}/${meta.tumor_id}.sig.allocation.tsv" : ''
    def cuppa_csv_arg = cuppa_dir ? "-cuppa_result_csv ${cuppa_dir}/${meta.tumor_id}.cup.data.csv" : ''
    def cuppa_summary_arg = cuppa_dir ? "-cuppa_summary_plot ${cuppa_dir}/${meta.tumor_id}.cup.report.summary.png" : ''

    def normal_id_arg = meta.containsKey('normal_wgs_id') ? "-reference_sample_id ${meta.normal_wgs_id}" : ''
    def normal_metrics_arg = bam_metrics_germline ? "-ref_sample_wgs_metrics_file ${bam_metrics_germline}" : ''
    def normal_flagstat_arg = flagstat_germline ? "-ref_sample_flagstat_file ${flagstat_germline}" : ''
    def normal_sage_somatic_bqr_arg = sage_germline_bqr ? "-sage_somatic_ref_sample_bqr_plot ${sage_germline_bqr}" : ''
    def normal_sage_coverage_arg = sage_germline_coverage ? "-sage_germline_gene_coverage_tsv ${sage_germline_coverage}" : ''
    def normal_linx_arg = linx_germline_anno_dir ? "-linx_germline_data_directory ${linx_germline_anno_dir}" : ''

    def rna_id_arg = meta.containsKey('tumor_wts_id') ? "-rna_sample_id ${meta.tumor_wts_id}" : ''
    def isofox_summary_csv_arg = isofox_dir ? "-isofox_summary_csv ${isofox_dir}/${meta.tumor_wts_id}.isf.summary.csv" : ''
    def isofox_gene_csv_arg = isofox_dir ? "-isofox_gene_data_csv ${isofox_dir}/${meta.tumor_wts_id}.isf.gene_data.csv" : ''
    def isofox_fusion_csv_arg = isofox_dir ? "-isofox_fusion_csv ${isofox_dir}/${meta.tumor_wts_id}.isf.fusions.csv" : ''
    def isofox_alt_sj_csv_arg = isofox_dir ? "-isofox_alt_splice_junction_csv ${isofox_dir}/${meta.tumor_wts_id}.isf.alt_splice_junc.csv" : ''

    def isofox_gene_distribution_arg = isofox_gene_distribution ? "-isofox_gene_distribution_csv ${isofox_gene_distribution}" : ''
    def isofox_alt_sj_arg = isofox_alt_sj ? "-isofox_alt_sj_cohort_csv ${isofox_alt_sj}" : ''

    """
    echo "${pipeline_version_str}" > pipeline_version.txt

    # When WTS data is present, ORANGE expects the somatic SAGE VCF to have appended WTS data; CS indicates this should
    # occur after PURPLE. Since ORANGE only collects the somatic SAGE VCF from the PURPLE output directory, we must
    # prepare accordingly
    purple_dir_local=${purple_dir}
    if [[ -n "${rna_id_arg}" ]]; then
        purple_dir_local=purple__prepared;
        mkdir -p \${purple_dir_local}/;
        find -L ${purple_dir} -maxdepth 1 -exec ln -fs ../{} \${purple_dir_local}/ \\;
        ln -sf ../${smlv_somatic_vcf} \${purple_dir_local}/${meta.tumor_id}.purple.somatic.vcf.gz;
        ln -sf ../${smlv_germline_vcf} \${purple_dir_local}/${meta.tumor_id}.purple.germline.vcf.gz;
    fi

    # NOTE(SW): '--add-opens java.base/java.time=ALL-UNNAMED' resolves issue writing JSON, see:
    # https://stackoverflow.com/questions/70412805/what-does-this-error-mean-java-lang-reflect-inaccessibleobjectexception-unable/70878195#70878195

    # NOTE(SW): DOID label: 162 [cancer]; Hartwig cohort group: unknown

    java \\
        --add-opens java.base/java.time=ALL-UNNAMED \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            \\
            -experiment_date \$(date +%y%m%d) \\
            -pipeline_version_file pipeline_version.txt \\
            \\
            -tumor_sample_id ${meta.tumor_id} \\
            -primary_tumor_doids 162 \\
            -tumor_sample_wgs_metrics_file ${bam_metrics_somatic} \\
            -tumor_sample_flagstat_file ${flagstat_somatic} \\
            -sage_somatic_tumor_sample_bqr_plot ${sage_somatic_bqr} \\
            -purple_data_directory \${purple_dir_local} \\
            -purple_plot_directory \${purple_dir_local}/plot/ \\
            -linx_somatic_data_directory ${linx_somatic_anno_dir} \\
            -linx_plot_directory ${linx_somatic_plot_dir} \\
            -lilac_result_csv ${lilac_dir}/${meta.tumor_id}.lilac.csv \\
            -lilac_qc_csv ${lilac_dir}/${meta.tumor_id}.lilac.qc.csv \\
            ${virusinterpreter_arg} \\
            ${chord_arg} \\
            ${sigs_arg} \\
            ${cuppa_csv_arg} \\
            ${cuppa_summary_arg} \\
            \\
            ${normal_id_arg} \\
            ${normal_metrics_arg} \\
            ${normal_flagstat_arg} \\
            ${normal_sage_somatic_bqr_arg} \\
            ${normal_sage_coverage_arg} \\
            ${normal_linx_arg} \\
            \\
            ${rna_id_arg} \\
            ${isofox_summary_csv_arg} \\
            ${isofox_gene_csv_arg} \\
            ${isofox_fusion_csv_arg} \\
            ${isofox_alt_sj_csv_arg} \\
            \\
            -ref_genome_version ${genome_ver} \\
            -doid_json ${disease_ontology} \\
            -cohort_mapping_tsv ${cohort_mapping} \\
            -cohort_percentiles_tsv ${cohort_percentiles} \\
            -known_fusion_file ${known_fusion_data} \\
            -driver_gene_panel_tsv ${driver_gene_panel} \\
            -ensembl_data_directory ${ensembl_data_resources} \\
            ${isofox_gene_distribution_arg} \\
            ${isofox_alt_sj_arg} \\
            -output_dir ./output/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orange: \$(java -jar ${task.ext.jarPath} | head -n1 | sed 's/.*ORANGE v//')
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
