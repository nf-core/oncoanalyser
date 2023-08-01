//
// CHORD predicts HR status for tumor samples
//

import Constants
import Utils

include { CHORD } from '../../modules/local/chord/main'

workflow CHORD_PREDICTION {
    take:
        // Sample data
        ch_inputs      // channel: [mandatory] [ meta ]
        ch_purple      // channel: [mandatory] [ meta, purple_dir ]

        // Reference data
        genome_version // channel: [mandatory] genome version

        // Params
        run_config     // channel: [mandatory] run configuration

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [ meta, purple_dir ]
        if (run_config.stages.purple) {
          ch_chord_inputs_source = ch_purple
        } else {
          ch_chord_inputs_source = WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PURPLE_DIR)
        }

        // Create inputs and create process-specific meta
        // channel: [ meta, smlv_vcf, sv_vcf ]
        ch_chord_inputs = ch_chord_inputs_source
            .map { meta, purple_dir ->
                def tumor_id = Utils.getTumorWgsSampleName(meta)
                def smlv_vcf = file(purple_dir).resolve("${tumor_id}.purple.somatic.vcf.gz")
                def sv_vcf = file(purple_dir).resolve("${tumor_id}.purple.sv.vcf.gz")

                // Require both SV and smlv VCF from the PURPLE directory
                if (!smlv_vcf.exists() || !sv_vcf.exists()) {
                    return Constants.PLACEHOLDER_META
                }

                def meta_chord = [key: meta.id, id: tumor_id]
                return [meta_chord, smlv_vcf, sv_vcf]
            }
            .filter { it != Constants.PLACEHOLDER_META }

        // Run process
        CHORD(
          ch_chord_inputs,
          genome_version,
        )

        // Set outputs, restoring original meta
        ch_outputs = WorkflowOncoanalyser.restoreMeta(CHORD.out.chord_dir, ch_inputs)
        ch_versions = ch_versions.mix(CHORD.out.versions)

    emit:
        chord_dir = ch_outputs // channel: [ meta, chord_dir ]

        versions = ch_versions // channel: [ versions.yml ]
}
