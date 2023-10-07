//
// Construct PURPLE input channel
//

import Constants

workflow CHANNEL_INPUTS_PURPLE {
    take:
        // Sample data
        ch_inputs                // channel: [mandatory] [ meta ]
        ch_amber                 // channel: [mandatory] [ meta, amber_dir ]
        ch_cobalt                // channel: [mandatory] [ meta, cobalt_dir ]
        ch_smlv_somatic          // channel: [optional]  [ meta, pave_vcf ]
        ch_smlv_germline         // channel: [optional]  [ meta, pave_vcf ]
        ch_sv_somatic            // channel: [optional]  [ meta, gripss_vcf, gripss_tbi ]
        ch_sv_germline           // channel: [optional]  [ meta, gripss_vcf, gripss_tbi ]
        ch_sv_somatic_unfiltered // channel: [optional]  [ meta, gripss_vcf, gripss_tbi ]

    main:

        ch_inputs_combined = WorkflowOncoanalyser.groupByMeta(
            ch_amber,
            ch_cobalt,
            ch_sv_somatic,
            ch_sv_somatic_unfiltered,
            ch_sv_germline,
            ch_smlv_somatic,
            ch_smlv_germline,
        )
            .map { d ->

                // Select sources here

            }



        /*

        // operating on basis of a runnable check
        //   * some logic to check whether required inputs are present
        //     * moving to strongly adhere, we could emit all non-runnable placeholders and merge/group here
        //     * then robust logic to check inputs to determine runnability
        //
        //
        //
        //       * replace non-empty input with samplesheet inputs
        //
        //         * DECIDE WHEN THIS IS DONE: TWO CASES:
        //              * stage run and stage NOT run
        //
        //
        //
        //       * then evaluate
        //
        //   * stages not run also need to be handled right at merge
        //     * this is where the run_config.stages.gripss (etc) are checked
        //     * initially intended to be handled by samplesheet retrieval, though now required earlier for merging
        //       * could an else block be added to call of the process subworkflow to generate these there?



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

        // NOTE(SW): Hartwig indicated unfiltered SVs should not be used for recovery in panel or tumor-only mode
        if (params.targeted === true || run_config.type == Constants.RunType.TUMOR_ONLY) {
            ch_sv_somatic_unfiltered_source = ch_inputs.map { meta -> [meta, [], []] }
        }

        // Combine into single channel
        ch_purple_inputs_source = WorkflowOncoanalyser.groupByMeta(
            ch_amber_source,
            ch_cobalt_source,
            ch_sv_somatic_source,
            ch_sv_somatic_unfiltered_source,
            ch_sv_germline_source,
            ch_smlv_somatic_source,
            ch_smlv_germline_source,
        )

        */

        ch_purple_inputs_source = Channel.empty()

    emit:
        ch_purple_inputs_source // channel: [ meta, amber_dir, cobalt_dir, sv_somatic_vcf, sv_somatic_tbi, sv_somatic_unfiltered_vcf, sv_somatic_unfiltered_tbi, sv_germline_vcf, sv_germline_tbi, smlv_somatic_vcf, smlv_germline_vcf ]
}

def getGripssSampleSheetInput(ch_inputs, input_type) {
    WorkflowOncoanalyser.getInput(ch_inputs, input_type, type: 'optional')
        .map { meta, vcf ->
            def tbi = vcf == [] ? [] : "${vcf}.tbi"
            return [meta, vcf, tbi]
        }
}
