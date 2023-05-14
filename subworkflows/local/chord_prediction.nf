//
// XXX
//
import Constants
import Utils

include { CHORD } from '../../modules/local/chord/main'

workflow CHORD_PREDICTION {
    take:
        // Sample data
        ch_inputs
        ch_purple

        // Reference data
        ref_data_genome_version

        // Parameters
        run

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [val(meta), purple_dir]
        if (run.purple) {
          ch_chord_inputs_source = ch_purple
        } else {
          ch_chord_inputs_source = WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PURPLE_DIR)
        }

        // Create inputs and create process-specific meta
        // channel: [val(meta), smlv_vcf, sv_vcf]
        ch_chord_inputs = ch_chord_inputs_source
            .map { meta, purple_dir ->
                def tumor_id = Utils.getTumorWgsSampleName(meta)
                def smlv_vcf = file(purple_dir).resolve("${tumor_id}.purple.somatic.vcf.gz")
                def sv_vcf = file(purple_dir).resolve("${tumor_id}.purple.sv.vcf.gz")

                // Require both SV and smlv VCF from the PURPLE directory
                if (!smlv_vcf.exists() || !sv_vcf.exists()) {
                    return Constants.META_PLACEHOLDER
                }

                def meta_chord = [key: meta.id, id: meta.id]
                return [meta_chord, smlv_vcf, sv_vcf]
            }
            .filter { it != Constants.META_PLACEHOLDER }

        // Run process
        CHORD(
          ch_chord_inputs,
          ref_data_genome_version,
        )

        // Set outputs, restoring original meta
        ch_outputs = WorkflowOncoanalyser.restoreMeta(CHORD.out.prediction, ch_inputs)
        ch_versions = ch_versions.mix(CHORD.out.versions)

    emit:
        prediction = ch_outputs // channel: [val(meta), prediction]

        versions = ch_versions  // channel: [versions.yml]
}

