//
// PEACH infers germline haplotypes and reports relevant pharmacogenomics
//

include { PEACH } from '../../../modules/local/peach/main'

workflow PEACH_CALLING {
    take:
    // Sample data
    ch_inputs                 // channel: [mandatory] [ meta ]
    ch_purple_dir             // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    peach_haplotypes          // channel: [mandatory] /path/to/peach_haplotypes
    peach_haplotype_functions // channel: [mandatory] /path/to/peach_haplotype_functions
    peach_drug_info           // channel: [mandatory] /path/to/peach_drug_info

    main:
    // Select input sources then sort
    // channel: runnable: [ meta, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_purple_dir
        .map { meta, purple_dir ->
            return [
                meta,
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
            ]
        }
        .branch { meta, purple_dir ->

            def has_normal = Utils.hasNormalDna(meta)
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.PEACH_DIR)

            def tumor_id = Utils.getTumorDnaSampleName(meta)
            def has_smlv_vcf = purple_dir ? purple_dir.resolve("${tumor_id}.purple.germline.vcf.gz").exists() : false

            runnable: has_smlv_vcf && has_normal && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_peach, purple_germline_smlv_vcf ]
    ch_peach_inputs = ch_inputs_sorted.runnable
        .map { meta, purple_dir ->

            def meta_peach = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getNormalDnaSampleName(meta),
            ]

            def purple_germline_smlv_vcf = purple_dir.resolve("${Utils.getTumorDnaSampleName(meta)}.purple.germline.vcf.gz")

            return [meta_peach, purple_germline_smlv_vcf]
        }

    // Run process
    PEACH(
        ch_peach_inputs,
        peach_haplotypes,
        peach_haplotype_functions,
        peach_drug_info,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, peach_dir ]
    ch_outputs = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('peach_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    peach_dir = ch_outputs // channel: [ meta, peach_dir ]
}
