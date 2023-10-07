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
        ch_smlv_somatic          // channel: [mandatory] [ meta, pave_vcf ]
        ch_smlv_germline         // channel: [mandatory] [ meta, pave_vcf ]
        ch_sv_somatic            // channel: [mandatory] [ meta, gripss_vcf, gripss_tbi ]
        ch_sv_germline           // channel: [mandatory] [ meta, gripss_vcf, gripss_tbi ]
        ch_sv_somatic_unfiltered // channel: [mandatory] [ meta, gripss_vcf, gripss_tbi ]

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

        ch_inputs_selected = ch_inputs_combined
            .map { d ->

                def meta = d[0]

                // NOTE(SW): avoiding further complexity with loops etc

                def inputs = [
                    Utils.selectCurrentOrExisting(d[1], meta, Constants.INPUT.AMBER_DIR),
                    Utils.selectCurrentOrExisting(d[2], meta, Constants.INPUT.COBALT_DIR),
                    Utils.selectCurrentOrExisting(d[3], meta, Constants.INPUT.GRIPSS_VCF_TUMOR),
                    Utils.selectCurrentOrExisting(d[4], meta, Constants.INPUT.GRIPSS_VCF_TUMOR_TBI),
                    Utils.selectCurrentOrExisting(d[5], meta, Constants.INPUT.GRIPSS_UNFILTERED_VCF_TUMOR),
                    Utils.selectCurrentOrExisting(d[6], meta, Constants.INPUT.GRIPSS_UNFILTERED_VCF_TUMOR_TBI),
                    Utils.selectCurrentOrExisting(d[7], meta, Constants.INPUT.GRIPSS_VCF_NORMAL),
                    Utils.selectCurrentOrExisting(d[8], meta, Constants.INPUT.GRIPSS_VCF_NORMAL_TBI),
                    Utils.selectCurrentOrExisting(d[9], meta, Constants.INPUT.PAVE_VCF_TUMOR),
                    Utils.selectCurrentOrExisting(d[10], meta, Constants.INPUT.PAVE_VCF_NORMAL),
                ]

                return [meta, *inputs]

            }

        ch_inputs_sorted = ch_inputs_selected
            .branch { d ->
                def meta = d[0]
                def amber_dir = d[1]
                def cobalt_dir = d[2]

                def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.PURPLE_DIR)

                runnable: amber_dir && cobalt_dir && !has_existing
                skip: true
                    return meta
            }

    emit:
        ch_purple_inputs_source // channel: [ meta, amber_dir, cobalt_dir, sv_somatic_vcf, sv_somatic_tbi, sv_somatic_unfiltered_vcf, sv_somatic_unfiltered_tbi, sv_germline_vcf, sv_germline_tbi, smlv_somatic_vcf, smlv_germline_vcf ]
}
