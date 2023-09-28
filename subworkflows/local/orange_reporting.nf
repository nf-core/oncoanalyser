//
// ORANGE collates outputs of hmftools into a static PDF report
//

import Constants
import Utils

include { ORANGE } from '../../modules/local/orange/main'

include { FLAGSTAT_METRICS } from './flagstat_metrics'

workflow ORANGE_REPORTING {
    take:
        // Sample data
        ch_inputs                   // channel: [mandatory] [ meta ]
        ch_bamtools_somatic         // channel: [mandatory] [ meta, metrics ]
        ch_bamtools_germline        // channel: [optional]  [ meta, metrics ]
        ch_sage_somatic             // channel: [mandatory] [ meta, sage_dir ]
        ch_sage_germline            // channel: [optional]  [ meta, sage_dir ]
        ch_sage_somatic_append      // channel: [optional]  [ meta, sage_append_vcf ]
        ch_sage_germline_append     // channel: [optional]  [ meta, sage_append_vcf ]
        ch_purple                   // channel: [mandatory] [ meta, purple_dir ]
        ch_linx_somatic_annotation  // channel: [mandatory] [ meta, linx_annotation_dir ]
        ch_linx_somatic_plot        // channel: [mandatory] [ meta, linx_visualiser_dir ]
        ch_linx_germline_annotation // channel: [optional]  [ meta, linx_annotation_dir ]
        ch_virusinterpreter         // channel: [optional]  [ meta, virusinterpreter_dir ]
        ch_chord                    // channel: [optional]  [ meta, chord_dir ]
        ch_sigs                     // channel: [optional]  [ meta, sigs_dir ]
        ch_lilac                    // channel: [mandatory] [ meta, lilac_dir ]
        ch_cuppa                    // channel: [optional]  [ meta, cuppa_dir ]
        ch_isofox                   // channel: [optional]  [ meta, isofox_dir ]

        // Reference data
        genome_version              // channel: [mandatory] genome version
        disease_ontology            // channel: [mandatory] /path/to/disease_ontology
        cohort_mapping              // channel: [mandatory] /path/to/cohort_mapping
        cohort_percentiles          // channel: [mandatory] /path/to/cohort_percentiles
        known_fusion_data           // channel: [mandatory] /path/to/known_fusion_data
        driver_gene_panel           // channel: [mandatory] /path/to/driver_gene_panel
        ensembl_data_resources      // channel: [mandatory] /path/to/ensembl_data_resources/
        isofox_alt_sj               // channel: [optional]  /path/to/isofox_alt_sj
        isofox_gene_distribution    // channel: [optional]  /path/to/isofox_gene_distribution

        // Params
        run_config                  // channel: [mandatory] run configuration

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        //
        // SUBWORKFLOW: Run SAMtools flagstat to generate stats required for ORANGE
        //
        // channel: [ meta, metrics ]
        ch_flagstat_somatic_out = Channel.empty()
        ch_flagstat_germline_out = Channel.empty()
        if (run_config.stages.flagstat) {

            FLAGSTAT_METRICS(
                ch_inputs,
                run_config,
            )

            ch_versions = ch_versions.mix(FLAGSTAT_METRICS.out.versions)
            ch_flagstat_somatic_out = ch_flagstat_somatic_out.mix(FLAGSTAT_METRICS.out.somatic)
            ch_flagstat_germline_out = ch_flagstat_germline_out.mix(FLAGSTAT_METRICS.out.germline)
        }

        //
        // MODULE: Run ORANGE
        //
        // Create placeholders for tumor-only
        if (run_config.type == Constants.RunType.TUMOR_ONLY) {
            ch_flagstat_germline_out = ch_inputs.map { meta -> [meta, []] }
            ch_bamtools_germline = ch_inputs.map { meta -> [meta, []] }
            ch_sage_germline = ch_inputs.map { meta -> [meta, []] }
            ch_sage_germline_append = ch_inputs.map { meta -> [meta, []] }
            ch_linx_germline_annotation = ch_inputs.map { meta -> [meta, []] }
        }

        // Do not pass WTS reference files when not required
        // NOTE(SW): ORANGE v2.6.0 will crash if WTS reference files but not WTS sample files provided
        if (run_config.mode != Constants.RunMode.WTS && run_config.mode != Constants.RunMode.WGTS) {
            isofox_alt_sj = []
            isofox_gene_distribution = []
        }

        // Get PURPLE input source for processing
        ch_orange_inputs_purple_dir = run_config.stages.purple ? ch_purple : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PURPLE_DIR)

        // Get input smlv somatic VCF from either PURPLE or SAGE append
        // channel: [ meta, sage_somatic_vcf, sage_germline_vcf ]
        if (run_config.mode == Constants.RunMode.WGS || run_config.mode == Constants.RunMode.PANEL) {

            ch_orange_inputs_smlv_vcfs = ch_orange_inputs_purple_dir
                .map { meta, purple_dir ->
                    def tumor_id = Utils.getTumorSampleName(meta, run_config.mode)

                    def smlv_somatic_vcf = []
                    def smlv_germline_vcf = []

                    def smlv_somatic_vcf_path = file(purple_dir).resolve("${tumor_id}.purple.somatic.vcf.gz")
                    def smlv_germline_vcf_path = file(purple_dir).resolve("${tumor_id}.purple.germline.vcf.gz")

                    if (smlv_somatic_vcf_path.exists()) {
                        smlv_somatic_vcf = smlv_somatic_vcf_path
                    }

                    // NOTE(SW): can only evaluate to true with WGS tumor/normal
                    if (smlv_germline_vcf_path.exists()) {
                        smlv_germline_vcf = smlv_germline_vcf_path
                    }

                    return [meta, smlv_somatic_vcf, smlv_germline_vcf]
                }

        } else if (run_config.mode == Constants.RunMode.WGTS) {

            ch_orange_inputs_smlv_vcfs = WorkflowOncoanalyser.groupByMeta(
                ch_sage_somatic_append,
                ch_sage_germline_append,
            )

        } else {
            ch_orange_inputs_smlv_vcfs = ch_inputs.map { meta -> [meta, [], []] }
        }

        // Select input source
        ch_orange_inputs_source = WorkflowOncoanalyser.groupByMeta(
            run_config.stages.bamtools ? ch_bamtools_somatic : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.BAMTOOLS_TUMOR),
            run_config.stages.bamtools ? ch_bamtools_germline : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.BAMTOOLS_NORMAL, type: 'optional'),
            run_config.stages.flagstat ? ch_flagstat_somatic_out : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.FLAGSTAT_TUMOR),
            run_config.stages.flagstat ? ch_flagstat_germline_out : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.FLAGSTAT_NORMAL, type: 'optional'),
            run_config.stages.sage ? ch_sage_somatic : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.SAGE_DIR_TUMOR),
            run_config.stages.sage ? ch_sage_germline : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.SAGE_DIR_NORMAL, type: 'optional'),
            ch_orange_inputs_purple_dir,
            ch_orange_inputs_smlv_vcfs,
            run_config.stages.linx ? ch_linx_somatic_annotation : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.LINX_ANNO_DIR_TUMOR),
            run_config.stages.linx ? ch_linx_somatic_plot : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.LINX_PLOT_DIR_TUMOR),
            run_config.stages.linx ? ch_linx_germline_annotation : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.LINX_ANNO_DIR_NORMAL, type: 'optional'),
            run_config.stages.virusinterpreter ? ch_virusinterpreter : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.VIRUSINTERPRETER_DIR, type: 'optional'),
            run_config.stages.chord ? ch_chord : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.CHORD_DIR, type: 'optional'),
            run_config.stages.sigs ? ch_sigs : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.SIGS_DIR, type: 'optional'),
            run_config.stages.lilac ? ch_lilac : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.LILAC_DIR),
            run_config.stages.cuppa ? ch_cuppa : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.CUPPA_DIR, type: 'optional'),
            run_config.stages.isofox ? ch_isofox : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.ISOFOX_DIR, type: 'optional'),
            flatten_mode: 'nonrecursive',
        )

        ch_orange_inputs = ch_orange_inputs_source
            .map {
                def meta = it[0]

                // NOTE(SW): these attributes are optional
                def normal_wgs_id = meta.getAt(['sample_name', Constants.SampleType.NORMAL, Constants.SequenceType.WGS])
                def tumor_wts_id = meta.getAt(['sample_name', Constants.SampleType.TUMOR, Constants.SequenceType.WTS])

                def meta_orange = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: Utils.getTumorSampleName(meta, run_config.mode),
                ]

                // Add optional identifiers to meta
                if (normal_wgs_id) {
                    meta_orange.normal_wgs_id = normal_wgs_id
                }
                if (tumor_wts_id) {
                    meta_orange.tumor_wts_id = tumor_wts_id
                }

                return [meta_orange, *it[1..-1]]
            }

        // Run process
        ORANGE(
            ch_orange_inputs,
            genome_version,
            disease_ontology,
            cohort_mapping,
            cohort_percentiles,
            known_fusion_data,
            driver_gene_panel,
            ensembl_data_resources,
            isofox_alt_sj,
            isofox_gene_distribution,
            "5.33 [oncoanalyser]",
        )

        // Set outputs
        ch_versions = ch_versions.mix(ORANGE.out.versions)

    emit:
        versions  = ch_versions // channel: [ versions.yml ]
}
