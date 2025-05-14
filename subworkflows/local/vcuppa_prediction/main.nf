//
// vCUPPA predicts tissue of origin from molecular profiles for targeted data
//

import Constants
import Utils

include { VCUPPA } from '../../../modules/local/vcuppa/main'

workflow VCUPPA_PREDICTION {
    take:
    // Sample data
    ch_inputs            // channel: [mandatory] [ meta ]
    ch_purple            // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_version       // channel: [mandatory] genome version
    vcuppa_model         // channel: [mandatory] /path/to/vcuppa_model
    vcuppa_features_list // channel: [mandatory] /path/to/vcuppa_features_list

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
            runnable: purple_dir
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_vcuppa, purple_dir ]
    ch_vcuppa_inputs = ch_inputs_sorted.runnable
        .map{ meta, purple_dir ->

            def meta_vcuppa = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Utils.getTumorDnaSampleName(meta),
            ]

            return [meta_vcuppa, purple_dir]
        }

    // Run process
    VCUPPA(
        ch_vcuppa_inputs,
        genome_version,
        vcuppa_model,
        vcuppa_features_list,
    )

    ch_versions = ch_versions.mix(VCUPPA.out.versions)

    emit:
    versions  = ch_versions // channel: [ versions.yml ]
}
