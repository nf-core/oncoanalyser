//
// Sigs fits trinucleotide signature definitions with sample SNV counts
//

import Constants
import Utils

include { SIGS } from '../../../modules/local/sigs/main'

workflow SIGS_FITTING {
    take:
    // Sample data
    ch_inputs       // channel: [mandatory] [ meta ]
    ch_purple       // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    sigs_signatures // channel: [mandatory] /path/to/sigs_signatures

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources
    // channel: [ meta, purple_dir ]
    ch_inputs_selected = ch_purple
        .map { meta, purple_dir ->
            return [meta, Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR)]
        }

    // Sort inputs
    // channel: runnable: [ meta, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_inputs_selected
        .branch { meta, purple_dir ->

            def has_dna = Utils.hasTumorDna(meta)

            def tumor_id
            def has_smlv_vcf
            if (has_dna) {
                tumor_id = Utils.getTumorDnaSampleName(meta)
                has_smlv_vcf = purple_dir ? file(purple_dir).resolve("${tumor_id}.purple.somatic.vcf.gz") : []
            }

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SIGS_DIR)

            runnable: has_dna && purple_dir && has_smlv_vcf && !has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_sigs, smlv_vcf ]
    ch_sigs_inputs = ch_inputs_sorted.runnable
        .map { meta, purple_dir ->

            def tumor_id = Utils.getTumorDnaSampleName(meta)

            def meta_sigs = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: tumor_id,
            ]

            def smlv_vcf = file(purple_dir).resolve("${tumor_id}.purple.somatic.vcf.gz")

            return [meta_sigs, smlv_vcf]
        }

    // Run process
    SIGS(
        ch_sigs_inputs,
        sigs_signatures,
    )

    ch_versions = ch_versions.mix(SIGS.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, sigs_dir ]
    ch_outputs = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(SIGS.out.sigs_dir, ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    sigs_dir = ch_outputs  // channel: [ meta, sigs_dir ]

    versions = ch_versions // channel: [ versions.yml ]
}
