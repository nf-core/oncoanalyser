//
// COBALT normalisation prepares the panel-specific target region normalisation resource
//

import Constants
import Utils

include { COBALT_PANEL_NORMALISATION } from '../../../modules/local/cobalt/panel_normalisation/main'

workflow COBALT_NORMALISATION {
    take:
    // Sample data
    ch_amber          // channel: [mandatory] [ meta, amber_dir ]
    ch_cobalt         // channel: [mandatory] [ meta, cobalt_dir ]

    // Reference data
    genome_version    // channel: [mandatory] genome version
    gc_profile        // channel: [mandatory] /path/to/gc_profile
    target_region_bed // channel: [mandatory] /path/to/target_region_bed

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Create process input channel
    // channel: [ [amber_dir, ...], [cobalt_dir, ...] ]
    ch_cobalt_inputs = WorkflowOncoanalyser.groupByMeta(
        ch_amber,
        ch_cobalt,
    )
        .map { meta, amber_dir, cobalt_dir ->
            return [
                Utils.selectCurrentOrExisting(amber_dir, meta, Constants.INPUT.AMBER_DIR),
                Utils.selectCurrentOrExisting(cobalt_dir, meta, Constants.INPUT.COBALT_DIR),
            ]
        }
        .collect(flat: false)
        .map { d -> d.transpose() }


    // Run process
    COBALT_PANEL_NORMALISATION(
        ch_cobalt_inputs,
        genome_version,
        gc_profile,
        target_region_bed,
    )

    ch_versions = ch_versions.mix(COBALT_PANEL_NORMALISATION.out.versions)

    emit:
    versions  = ch_versions // channel: [ versions.yml ]
}
