process ORANGE {
    tag "${meta.id}"
    label 'process_single'

    container 'docker.io/scwatts/orange:2.3--0'

    input:
    tuple val(meta), path(bam_metrics_somatic), path(bam_metrics_germline), path(flagstat_somatic), path(flagstat_germline), path(sage_somatic_bqr), path(sage_germline_bqr), path(sage_germline_coverage), path(purple_dir), path(linx_somatic_anno_dir), path(linx_somatic_plot_dir), path(linx_germline_anno_dir), path(virusinterpreter), path(chord_prediction), path(sigs_dir), path(lilac_dir), path(cuppa_dir), path(isofox_dir)
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
    tuple val(meta), path('*.orange.pdf') , emit: pdf
    tuple val(meta), path('*.orange.json'), emit: json
    path 'versions.yml'                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def pipeline_version_str = pipeline_version ?: 'not specified'

    def virusinterpreter_arg = virusinterpreter ? "-annotated_virus_tsv ${virusinterpreter}" : ''
    def chord_arg = chord_prediction ? "-chord_prediction_txt ${chord_prediction}" : ''
    def sigs_arg = sigs_dir ? "-sigs_allocation_tsv ${sigs_dir}/${meta.tumor_wgs_id}.sig.allocation.tsv" : ''
    def cuppa_csv_arg = cuppa_dir ? "-cuppa_result_csv ${cuppa_dir}/${meta.tumor_wgs_id}.cup.data.csv" : ''
    def cuppa_summary_arg = cuppa_dir ? "-cuppa_summary_plot ${cuppa_dir}/${meta.tumor_wgs_id}.cup.report.summary.png" : ''

    def normal_id_arg = meta.containsKey('normal_wgs_id') ? "-reference_sample_id ${meta.normal_wgs_id}" : ''
    def normal_metrics_arg = bam_metrics_germline ? "-ref_sample_wgs_metrics_file ${bam_metrics_germline}" : ''
    def normal_flagstat_arg = flagstat_germline ? "-ref_sample_flagstat_file ${flagstat_germline}" : ''
    def normal_sage_somatic_bqr_arg = sage_germline_bqr ? "-sage_somatic_ref_sample_bqr_plot ${sage_germline_bqr}" : ''
    def normal_sage_coverage_arg = sage_germline_coverage ? "-sage_germline_gene_coverage_tsv ${sage_germline_coverage}" : ''
    // NOTE(SW): LINX germline outputs are processed below to contain tumor id in the filename as required by ORANGE
    def normal_linx_arg = linx_germline_anno_dir ? "-linx_germline_data_directory normal_linx__prepared/" : ''

    def rna_id_arg = meta.containsKey('tumor_wts_id') ? "-rna_sample_id ${meta.tumor_wts_id}" : ''
    def isofox_summary_csv_arg = isofox_dir ? "-isofox_summary_csv ${isofox_dir}/${meta.tumor_wts_id}.isf.summary.csv" : ''
    def isofox_gene_csv_arg = isofox_dir ? "-isofox_gene_data_csv ${isofox_dir}/${meta.tumor_wts_id}.isf.gene_data.csv" : ''
    def isofox_fusion_csv_arg = isofox_dir ? "-isofox_fusion_csv ${isofox_dir}/${meta.tumor_wts_id}.isf.fusions.csv" : ''
    def isofox_alt_sj_csv_arg = isofox_dir ? "-isofox_alt_splice_junction_csv ${isofox_dir}/${meta.tumor_wts_id}.isf.alt_splice_junc.csv" : ''

    def isofox_gene_distribution_arg = isofox_gene_distribution ? "-isofox_gene_distribution_csv ${isofox_gene_distribution}" : ''
    def isofox_alt_sj_arg = isofox_alt_sj ? "-isofox_alt_sj_cohort_csv ${isofox_alt_sj}" : ''

    """
    echo "${pipeline_version_str}" > pipeline_version.txt

    if [[ -n "${normal_linx_arg}" ]]; then
        mkdir -p normal_linx__prepared/;
        for fp in ${linx_germline_anno_dir}/*; do
            fn=\${fp##*/};
            ln -s ../${linx_germline_anno_dir}/\${fn} normal_linx__prepared/\$(sed 's/${meta.normal_wgs_id}/${meta.tumor_wgs_id}/' <<< \${fn});
        done;
    fi

    # NOTE(SW): '--add-opens java.base/java.time=ALL-UNNAMED' resolves issue writing JSON, see:
    # https://stackoverflow.com/questions/70412805/what-does-this-error-mean-java-lang-reflect-inaccessibleobjectexception-unable/70878195#70878195

    java \\
        --add-opens java.base/java.time=ALL-UNNAMED \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -jar ${task.ext.jarPath} \\
            \\
            -experiment_date \$(date +%y%m%d) \\
            -pipeline_version_file pipeline_version.txt \\
            \\
            -tumor_sample_id ${meta.tumor_wgs_id} \\
            -primary_tumor_doids "" \\
            -tumor_sample_wgs_metrics_file ${bam_metrics_somatic} \\
            -tumor_sample_flagstat_file ${flagstat_somatic} \\
            -sage_somatic_tumor_sample_bqr_plot ${sage_somatic_bqr} \\
            -purple_data_directory ${purple_dir} \\
            -purple_plot_directory ${purple_dir}/plot/ \\
            -linx_somatic_data_directory ${linx_somatic_anno_dir} \\
            -linx_plot_directory ${linx_somatic_plot_dir} \\
            -lilac_result_csv ${lilac_dir}/${meta.tumor_wgs_id}.lilac.csv \\
            -lilac_qc_csv ${lilac_dir}/${meta.tumor_wgs_id}.lilac.qc.csv \\
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
            \\ #${rna_id_arg} \\
            \\ #${isofox_summary_csv_arg} \\
            \\ #${isofox_gene_csv_arg} \\
            \\ #${isofox_fusion_csv_arg} \\
            \\ #${isofox_alt_sj_csv_arg} \\
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
            -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orange: \$(java -jar ${task.ext.jarPath} | head -n1 | sed 's/.*ORANGE v//')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.tumor_wgs_id}.orange.json"
    touch "${meta.tumor_wgs_id}.orange.pdf"
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
