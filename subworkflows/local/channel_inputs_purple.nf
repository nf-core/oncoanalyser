//
// XXX
//
import Constants

workflow CHANNEL_INPUTS_PURPLE {
    take:
        // Sample data
        ch_inputs
        ch_amber
        ch_cobalt
        ch_smlv_somatic
        ch_smlv_germline
        ch_sv_somatic
        ch_sv_germline
        ch_sv_somatic_unfiltered

        // Params
        run_config

    main:
        // Set input sources
        ch_amber_source = run_config.stages.amber ? ch_amber : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.AMBER_DIR)
        ch_cobalt_source = run_config.stages.cobalt ? ch_cobalt : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.COBALT_DIR)

        if (run_config.stages.gripss) {

            ch_sv_somatic_source = ch_sv_somatic
            ch_sv_somatic_unfiltered_source = ch_sv_somatic_unfiltered
            ch_sv_germline_source = ch_sv_germline

        } else {

            ch_sv_somatic_source = getGripssSampleSheetInput(ch_inputs, Constants.INPUT.GRIPSS_VCF_TUMOR)
            ch_sv_somatic_unfiltered_source = getGripssSampleSheetInput(ch_inputs, Constants.INPUT.GRIPSS_UNFILTERED_VCF_TUMOR)
            ch_sv_germline_source = getGripssSampleSheetInput(ch_inputs, Constants.INPUT.GRIPSS_VCF_NORMAL)

        }

        ch_smlv_somatic_source = run_config.stages.pave ? ch_smlv_somatic : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PAVE_VCF_TUMOR, type: 'optional')
        ch_smlv_germline_source = run_config.stages.pave ? ch_smlv_germline : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PAVE_VCF_NORMAL, type: 'optional')

        // Adjust for run type and mode
        if (run_config.type == Constants.RunType.TUMOR_ONLY) {

            ch_sv_germline_source = ch_inputs.map { meta -> [meta, [], []] }
            ch_smlv_germline_source = ch_inputs.map { meta -> [meta, []] }

        }

        // Combine into single channel
        // channel: [val(meta), amber_dir, cobalt_dir, sv_tumor_vcf, sv_tumor_tbi, sv_tumor_unfiltered_vcf, sv_tumor_unfiltered_tbi, smlv_tumor_vcf, smlv_normal_vcf]
        ch_purple_inputs_source = WorkflowOncoanalyser.groupByMeta(
            ch_amber_source,
            ch_cobalt_source,
            ch_sv_somatic_source,
            ch_sv_somatic_unfiltered_source,
            ch_sv_germline_source,
            ch_smlv_somatic_source,
            ch_smlv_germline_source,
        )

    emit:
        ch_purple_inputs_source
}

def getGripssSampleSheetInput(ch_inputs, input_type) {
    WorkflowOncoanalyser.getInput(ch_inputs, input_type, type: 'optional')
        .map { meta, vcf ->
            def tbi = vcf == [] ? [] : "${vcf}.tbi"
            return [meta, vcf, tbi]
        }
}
