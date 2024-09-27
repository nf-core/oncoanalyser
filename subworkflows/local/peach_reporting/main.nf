//
// Pharmacogenomic Evaluator And Caller of Haplotypes (PEACH)
//

import Constants
import Utils

include { PEACH } from '../../../modules/local/peach/main'

workflow PEACH_REPORTING {
    take:
    ch_inputs               // channel: [mandatory] [ meta ]
    ch_sage_germline_vcf    // channel: [mandatory] [ meta, sage_germline_vcf, sage_somatic_tbi ]
    ch_haplotypes           // channel: [mandatory] [ meta, haplotypes_tsv ]
    ch_haplotype_functions  // channel: [mandatory] [ meta, haplotype_functions_tsv ]
    ch_drugs                // channel: [mandatory] [ meta, peach_drugs_tsv ]

    main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    //
    // MODULE: PAVE germline
    //
    // Select input sources and sort
    // channel: runnable: [ meta, sage_vcf, sage_tbi ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_sage_germline_vcf
        .map { meta, sage_vcf, sage_tbi ->
            return [
                meta,
                Utils.selectCurrentOrExisting(sage_vcf, meta, Constants.INPUT.SAGE_VCF_NORMAL),
                Utils.selectCurrentOrExisting(sage_tbi, meta, Constants.INPUT.SAGE_VCF_TBI_NORMAL),
            ]
        }
        .branch { meta, sage_vcf, sage_tbi ->

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.PEACH_DIR_NORMAL)

            runnable: Utils.hasTumorDna(meta) && Utils.hasNormalDna(meta) && sage_vcf && !has_existing
            skip: true
            return meta
        }

    // Create process input channel
    // channel: [ meta_pave, sage_vcf, sage_tbi ]
    ch_peach_inputs = ch_inputs_sorted.runnable
        .map { meta, sage_vcf, sage_tbi ->

            def meta_peach = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
            ]

            return [meta_peach, sage_vcf, sage_tbi]
        }

    // Run process
    PEACH(
        ch_peach_inputs,
        ch_haplotypes,
        ch_haplotype_functions,
        ch_drugs,
    )

    ch_versions = ch_versions.mix(PEACH.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, peach_dir ]
    ch_outputs = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(PEACH.out.peach_dir, ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    peach_dir = ch_outputs  // channel: [ meta, peach_dir ]

    versions = ch_versions  // channel: [ versions.yml ]
}
