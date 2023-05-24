//
// XXX
//
import Constants
import Utils

include { ORANGE } from '../../modules/local/orange/main'

include { CHANNEL_GROUP_INPUTS } from './channel_group_inputs'
include { FLAGSTAT_METRICS     } from './flagstat_metrics'

workflow ORANGE_REPORTING {
    take:
        // Sample data
        ch_inputs
        ch_bamtools_somatic
        ch_bamtools_germline
        ch_sage_somatic_tumor_bqr
        ch_sage_somatic_normal_bqr
        ch_sage_germline_coverage
        ch_sage_somatic_append
        ch_sage_germline_append
        ch_purple
        ch_linx_somatic_annotation
        ch_linx_somatic_plot
        ch_linx_germline_annotation
        ch_virusinterpreter
        ch_chord
        ch_sigs
        ch_lilac
        ch_cuppa
        ch_isofox

        // Reference data
        ref_data_genome_version
        ref_data_disease_ontology
        ref_data_cohort_mapping
        ref_data_cohort_percentiles
        ref_data_known_fusion_data
        ref_data_driver_gene_panel
        ref_data_ensembl_data_resources
        ref_data_isofox_alt_sj
        ref_data_isofox_gene_distribution

        // Params
        run

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Get input meta groups
        CHANNEL_GROUP_INPUTS(
            ch_inputs,
        )

        //
        // SUBWORKFLOW: Run SAMtools flagstat to generate stats required for ORANGE
        //
        // channel: [val(meta), metrics]
        ch_flagstat_somatic_out = Channel.empty()
        ch_flagstat_germline_out = Channel.empty()
        if (run.flagstat) {

            FLAGSTAT_METRICS(
                CHANNEL_GROUP_INPUTS.out.wgs_present,
            )

            ch_versions = ch_versions.mix(FLAGSTAT_METRICS.out.versions)
            ch_flagstat_somatic_out = ch_flagstat_somatic_out.mix(FLAGSTAT_METRICS.out.somatic)
            ch_flagstat_germline_out = ch_flagstat_germline_out.mix(FLAGSTAT_METRICS.out.germline)
        }

        //
        // MODULE: Run ORANGE
        //
        // Get PURPLE input source for processing
        ch_orange_inputs_purple_dir = run.purple ? ch_purple : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PURPLE_DIR)

        // Get input smlv somatic VCF from either PURPLE or SAGE append
        // NOTE(SW): a better test would be WTS or WGTS; this will be implemented as run mode (WGS/WTS/WGTS,
        // panel, MRD) along with input type (tumor-only, tumor/normal)
        // channel: [meta, smlv_somatic_vcf, smlv_germline_vcf]
        ch_orange_inputs_smlv_vcfs_append = WorkflowOncoanalyser.groupByMeta(
            ch_sage_somatic_append,
            ch_sage_germline_append,
        )
            // NOTE(SW): using .join for selection only
            .join(CHANNEL_GROUP_INPUTS.out.wts_present)

        // channel: [meta, smlv_somatic_vcf, smlv_germline_vcf]
        ch_orange_inputs_smlv_vcfs_purple_dir = CHANNEL_GROUP_INPUTS.out.wts_absent
            .join(ch_orange_inputs_purple_dir)
            .map { meta, purple_dir ->
                def tumor_id = Utils.getTumorWgsSampleName(meta)
                def smlv_somatic_vcf = file(purple_dir).resolve("${tumor_id}.purple.somatic.vcf.gz")
                def smlv_germline_vcf = file(purple_dir).resolve("${tumor_id}.purple.germline.vcf.gz")

                // Require smlv somatic VCF from the PURPLE directory
                if (!smlv_somatic_vcf.exists() || !smlv_germline_vcf.exists()) {
                    return Constants.META_PLACEHOLDER
                }

                return [meta, smlv_somatic_vcf, smlv_germline_vcf]
            }

        ch_orange_inputs_smlv_vcfs = Channel.empty()
            .mix(
                ch_orange_inputs_smlv_vcfs_append,
                ch_orange_inputs_smlv_vcfs_purple_dir,
            )

        // Select input source
        ch_orange_inputs_source = WorkflowOncoanalyser.groupByMeta(
            run.bamtools ? ch_bamtools_somatic : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.BAMTOOLS_TUMOR),
            run.bamtools ? ch_bamtools_germline : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.BAMTOOLS_NORMAL, type: 'optional'),
            run.flagstat ? ch_flagstat_somatic_out : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.FLAGSTAT_TUMOR),
            run.flagstat ? ch_flagstat_germline_out : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.FLAGSTAT_NORMAL, type: 'optional'),
            run.sage ? ch_sage_somatic_tumor_bqr : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.SAGE_BQR_TUMOR),
            run.sage ? ch_sage_somatic_normal_bqr : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.SAGE_BQR_NORMAL, type: 'optional'),
            run.sage ? ch_sage_germline_coverage : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.SAGE_COVERAGE, type: 'optional'),
            ch_orange_inputs_purple_dir,
            ch_orange_inputs_smlv_vcfs,
            run.linx ? ch_linx_somatic_annotation : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.LINX_ANNO_DIR_TUMOR),
            run.linx ? ch_linx_somatic_plot : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.LINX_PLOT_DIR_TUMOR),
            run.linx ? ch_linx_germline_annotation : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.LINX_PLOT_DIR_NORMAL, type: 'optional'),
            run.virusinterpreter ? ch_virusinterpreter : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.VIRUSINTERPRETER_TSV, type: 'optional'),
            run.chord ? ch_chord : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.CHORD_PREDICITION, type: 'optional'),
            run.sigs ? ch_sigs : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.SIGS_DIR, type: 'optional'),
            run.lilac ? ch_lilac : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.LILAC_DIR),
            run.cuppa ? ch_cuppa : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.CUPPA_DIR, type: 'optional'),
            run.isofox ? ch_isofox : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.ISOFOX_DIR, type: 'optional'),
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
                    tumor_wgs_id: Utils.getTumorWgsSampleName(meta),
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
            ref_data_genome_version,
            ref_data_disease_ontology,
            ref_data_cohort_mapping,
            ref_data_cohort_percentiles,
            ref_data_known_fusion_data,
            ref_data_driver_gene_panel,
            ref_data_ensembl_data_resources,
            ref_data_isofox_alt_sj,
            ref_data_isofox_gene_distribution,
            "5.32 [oncoanalyser]",
        )

        // Set outputs
        ch_versions = ch_versions.mix(ORANGE.out.versions)

    emit:
        versions  = ch_versions // channel: [versions.yml]
}
