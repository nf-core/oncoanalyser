//
// CHORD predicts HR status for tumor samples
//

import Constants
import Utils

include { CHORD } from '../../../modules/local/chord/main'

workflow CHORD_PREDICTION {
    take:
    // Sample data
    ch_inputs    // channel: [mandatory] [ meta ]
    ch_purple    // channel: [mandatory] [ meta, purple_dir ]
    genome_fasta // channel: [mandatory] /path/to/genome_fasta
    genome_fai   // channel: [mandatory] /path/to/genome_fai
    genome_dict  // channel: [mandatory] /path/to/genome_dict

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
            def has_sv_vcf

            if (has_dna) {
                tumor_id = Utils.getTumorDnaSampleName(meta)
                has_smlv_vcf = purple_dir ? file(purple_dir).resolve("${tumor_id}.purple.somatic.vcf.gz") : []
                has_sv_vcf = purple_dir ? file(purple_dir).resolve("${tumor_id}.purple.sv.vcf.gz") : []
            }

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.CHORD_DIR)

            runnable: has_dna && purple_dir && has_smlv_vcf && has_sv_vcf && !has_existing
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

            def smlv_vcf = file(purple_dir).resolve("${tumor_id}.purple.somatic.vcf.gz")
            def sv_vcf = file(purple_dir).resolve("${tumor_id}.purple.sv.vcf.gz")

            return [meta_chord, smlv_vcf, sv_vcf]
        }

    // Run process
    CHORD(
        ch_chord_inputs,
        genome_fasta,
        genome_fai,
        genome_dict,
    )

    ch_versions = ch_versions.mix(CHORD.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, chord_dir ]
    ch_outputs = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(CHORD.out.chord_dir, ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    chord_dir = ch_outputs  // channel: [ meta, chord_dir ]

    versions  = ch_versions // channel: [ versions.yml ]
}
