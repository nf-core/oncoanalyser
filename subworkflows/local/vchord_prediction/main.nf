//
// vCHORD predicts HR status for tumor samples from targeted data
//

import Constants
import Utils

include { VCHORD } from '../../../modules/local/vchord/main'

workflow VCHORD_PREDICTION {
    take:
    // Sample data
    ch_inputs    // channel: [mandatory] [ meta ]
    ch_purple    // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    vchord_model // channel: [mandatory] /path/to/vchord_model

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()


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
    // channel: [ meta_vchord, purple_dir ]
    ch_vchord_inputs = ch_inputs_sorted.runnable
        .map{ meta, purple_dir ->

            def meta_vchord = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Utils.getTumorDnaSampleName(meta),
            ]

            return [meta_vchord, purple_dir]
        }

    // Run process
    VCHORD(
        ch_vchord_inputs,
        vchord_model,
    )

    ch_versions = ch_versions.mix(VCHORD.out.versions)

    emit:
    versions  = ch_versions // channel: [ versions.yml ]
}
