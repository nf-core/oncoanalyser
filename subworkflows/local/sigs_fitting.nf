//
// Sigs fits trinucleotide signature definitions with sample SNV counts
//

import Constants
import Utils

include { SIGS } from '../../modules/local/sigs/main'

workflow SIGS_FITTING {
    take:
        // Sample data
        ch_inputs       // channel: [mandatory] [ meta ]
        ch_purple       // channel: [mandatory] [ meta, purple_dir ]

        // Reference data
        sigs_signatures // channel: [mandatory] /path/to/sigs_signatures

        // Params
        run_config      // channel: [mandatory] run configuration

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [ meta, purple_dir ]
        if (run_config.stages.purple) {
            ch_sigs_inputs_source = ch_purple
        } else {
            ch_sigs_inputs_source = WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PURPLE_DIR)
        }

        // Create inputs and create process-specific meta
        // channel: [ meta_sigs, smlv_vcf ]
        ch_sigs_inputs = ch_sigs_inputs_source
            .map { meta, purple_dir ->

                def tumor_id = Utils.getTumorDnaSampleName(meta)
                def smlv_vcf = file(purple_dir).resolve("${tumor_id}.purple.somatic.vcf.gz")

                // Require smlv VCF from the PURPLE directory
                if (!smlv_vcf.exists()) {
                    return Constants.PLACEHOLDER_META
                }

                def meta_sigs = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: tumor_id,
                ]
                return [meta_sigs, smlv_vcf]
            }
            .filter { it != Constants.PLACEHOLDER_META }

        SIGS(
          ch_sigs_inputs,
          sigs_signatures,
        )

        // Set outputs, restoring original meta
        ch_outputs = WorkflowOncoanalyser.restoreMeta(SIGS.out.sigs_dir, ch_inputs)
        ch_versions = ch_versions.mix(SIGS.out.versions)

    emit:
        sigs_dir = ch_outputs  // channel: [ meta, sigs_dir ]

        versions = ch_versions // channel: [ versions.yml ]
}
