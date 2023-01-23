// NOTE(SW): reference files for Isfox aren't currently available (-isofox_gene_distribution_csv, -isofox_alt_sj_cohort_csv)

process ORANGE {
    tag "${meta.id}"
    label 'process_single'

    container 'docker.io/scwatts/orange:1.10.2--0'

    input:
    tuple val(meta), path(bam_metrics_somatic), path(bam_metrics_germline), path(flagstat_somatic), path(flagstat_germline), path(chord_prediction), path(lilac_dir), path(sage_somatic_bqr), path(sage_germline_bqr), path(sage_germline_coverage), path(purple_dir), path(linx_somatic_anno_dir), path(linx_somatic_plot_dir), path(linx_germline_anno_dir), path(protect), path(peach_genotype), path(cuppa), path(cuppa_summary_plot), path(cuppa_feature_plot), path(virusinterpreter)
    val genome_ver
    path disease_ontology
    path known_fusion_data
    path driver_gene_panel
    path cohort_mapping
    path cohort_percentiles
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

    """
    echo "${pipeline_version_str}" > pipeline_version.txt

    # NOTE(SW): '--add-opens java.base/java.time=ALL-UNNAMED' resolves issue writing JSON, see:
    # https://stackoverflow.com/questions/70412805/what-does-this-error-mean-java-lang-reflect-inaccessibleobjectexception-unable/70878195#70878195

    java \\
        --add-opens java.base/java.time=ALL-UNNAMED \\
        -Xmx${task.memory.giga}g \\
        -jar ${task.ext.jarPath} \\
            \\
            -tumor_sample_id ${meta.tumor_id} \\
            -reference_sample_id ${meta.normal_id} \\
            \\
            -pipeline_version_file pipeline_version.txt \\
            -ref_genome_version ${genome_ver} \\
            -primary_tumor_doids "" \\
            -max_evidence_level C \\
            \\
            -tumor_sample_wgs_metrics_file ${bam_metrics_somatic} \\
            -ref_sample_wgs_metrics_file ${bam_metrics_germline} \\
            -tumor_sample_flagstat_file ${flagstat_somatic} \\
            -ref_sample_flagstat_file ${flagstat_germline} \\
            \\
            -chord_prediction_txt ${chord_prediction} \\
            -lilac_result_csv ${lilac_dir}/${meta.tumor_id}.lilac.csv \\
            -lilac_qc_csv ${lilac_dir}/${meta.tumor_id}.lilac.qc.csv \\
            \\
            -sage_somatic_tumor_sample_bqr_plot ${sage_somatic_bqr} \\
            -sage_somatic_ref_sample_bqr_plot ${sage_germline_bqr} \\
            -sage_germline_gene_coverage_tsv ${sage_germline_coverage} \\
            \\
            -purple_qc_file ${purple_dir}/${meta.tumor_id}.purple.qc \\
            -purple_purity_tsv ${purple_dir}/${meta.tumor_id}.purple.purity.tsv \\
            -purple_gene_copy_number_tsv ${purple_dir}/${meta.tumor_id}.purple.cnv.gene.tsv \\
            -purple_somatic_copy_number_tsv ${purple_dir}/${meta.tumor_id}.purple.cnv.somatic.tsv \\
            -purple_somatic_driver_catalog_tsv ${purple_dir}/${meta.tumor_id}.driver.catalog.somatic.tsv \\
            -purple_somatic_variant_vcf ${purple_dir}/${meta.tumor_id}.purple.somatic.vcf.gz \\
            -purple_germline_driver_catalog_tsv ${purple_dir}/${meta.tumor_id}.driver.catalog.germline.tsv \\
            -purple_germline_deletion_tsv ${purple_dir}/${meta.tumor_id}.purple.germline.deletion.tsv \\
            -purple_germline_variant_vcf ${purple_dir}/${meta.tumor_id}.purple.germline.vcf.gz \\
            -purple_plot_directory ${purple_dir}/plot/ \\
            \\
            -linx_breakend_tsv ${linx_somatic_anno_dir}/${meta.tumor_id}.linx.breakend.tsv \\
            -linx_structural_variant_tsv ${linx_somatic_anno_dir}/${meta.tumor_id}.linx.svs.tsv \\
            -linx_driver_tsv ${linx_somatic_anno_dir}/${meta.tumor_id}.linx.drivers.tsv \\
            -linx_driver_catalog_tsv ${linx_somatic_anno_dir}/${meta.tumor_id}.linx.driver.catalog.tsv \\
            -linx_fusion_tsv ${linx_somatic_anno_dir}/${meta.tumor_id}.linx.fusion.tsv \\
            -linx_germline_disruption_tsv ${linx_germline_anno_dir}/${meta.normal_id}.linx.germline.disruption.tsv \\
            -linx_plot_directory ${linx_somatic_plot_dir} \\
            \\
            -cuppa_result_csv ${cuppa} \\
            -cuppa_feature_plot ${cuppa_feature_plot} \\
            -cuppa_summary_plot ${cuppa_summary_plot} \\
            \\
            -peach_genotype_tsv ${peach_genotype} \\
            -protect_evidence_tsv ${protect} \\
            -annotated_virus_tsv ${virusinterpreter} \\
            \\
            -doid_json ${disease_ontology} \\
            -known_fusion_file ${known_fusion_data} \\
            -driver_gene_panel_tsv ${driver_gene_panel} \\
            -cohort_mapping_tsv ${cohort_mapping} \\
            -cohort_percentiles_tsv ${cohort_percentiles} \\
            -output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orange: \$(java -jar ${task.ext.jarPath} | head -n1 | sed 's/.*ORANGE v//')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.tumor_id}.orange.json"
    touch "${meta.tumor_id}.orange.pdf"
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
