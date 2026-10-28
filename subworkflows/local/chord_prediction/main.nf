//
// CHORD predicts HR status for tumor samples
//

include { CHORD } from '../../../modules/local/chord/main'

workflow CHORD_PREDICTION {
    take:
    // Sample data
    ch_inputs      // channel: [mandatory] [ meta ]
    ch_purple_dir  // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_fasta   // channel: [mandatory] /path/to/genome_fasta
    genome_fai     // channel: [mandatory] /path/to/genome_fai
    genome_dict    // channel: [mandatory] /path/to/genome_dict

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
            def has_sv_vcf = []
            if (has_tumor_normal_dna && purple_dir) {
                def tumor_id = Utils.getTumorDnaSampleName(meta)
                has_smlv_vcf = purple_dir.resolve("${tumor_id}.purple.somatic.vcf.gz").exists()
                has_sv_vcf = purple_dir.resolve("${tumor_id}.purple.sv.vcf.gz").exists()
            }

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.CHORD_DIR)

            runnable: has_tumor_normal_dna && purple_dir && has_smlv_vcf && has_sv_vcf && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_chord, smlv_vcf, sv_vcf ]
    ch_chord_inputs = ch_inputs_sorted.runnable
        .map { meta, purple_dir ->

            def tumor_id = Utils.getTumorDnaSampleName(meta)

            def meta_chord = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: tumor_id,
            ]

            def smlv_vcf = purple_dir.resolve("${tumor_id}.purple.somatic.vcf.gz")
            def sv_vcf = purple_dir.resolve("${tumor_id}.purple.sv.vcf.gz")

            return [meta_chord, smlv_vcf, sv_vcf]
        }

    // Run process
    CHORD(
        ch_chord_inputs,
        genome_fasta,
        genome_fai,
        genome_dict,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, chord_dir ]
    ch_outputs = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('chord_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    chord_dir = ch_outputs // channel: [ meta, chord_dir ]
}
