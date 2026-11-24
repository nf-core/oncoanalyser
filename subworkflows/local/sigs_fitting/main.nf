//
// Sigs fits trinucleotide signature definitions with sample SNV counts
//

include { SIGS } from '../../../modules/local/sigs/main'

workflow SIGS_FITTING {
    take:
    // Sample data
    ch_inputs       // channel: [mandatory] [ meta ]
    ch_purple_dir   // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    sigs_signatures // channel: [mandatory] /path/to/sigs_signatures

    main:
    // Select input sources then sort
    // channel: runnable: [ meta, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_purple_dir
        .map { meta, purple_dir ->
            return [meta, Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR)]
        }
        .branch { meta, purple_dir ->

            def has_tumor_normal_dna = Utils.hasTumorDna(meta) && Utils.hasNormalDna(meta)

            def has_smlv_vcf = []
            if (has_tumor_normal_dna && purple_dir) {
                def tumor_id = Utils.getTumorDnaSampleName(meta)
                has_smlv_vcf = purple_dir.resolve("${tumor_id}.purple.somatic.vcf.gz").exists()
            }

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SIGS_DIR)

            runnable: has_tumor_normal_dna && purple_dir && has_smlv_vcf && ! has_existing
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

            def smlv_vcf = purple_dir.resolve("${tumor_id}.purple.somatic.vcf.gz")

            return [meta_sigs, smlv_vcf]
        }

    // Run process
    SIGS(
        ch_sigs_inputs,
        sigs_signatures,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, sigs_dir ]
    ch_outputs = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('sigs_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    sigs_dir = ch_outputs // channel: [ meta, sigs_dir ]
}
