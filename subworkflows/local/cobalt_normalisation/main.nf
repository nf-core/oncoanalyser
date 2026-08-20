//
// COBALT normalisation prepares the panel-specific target region normalisation resource
//

include { COBALT_PANEL_NORMALISATION } from '../../../modules/local/cobalt/panel_normalisation/main'
include { groupByMeta             } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { joinMeta                } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { restoreMeta             } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { getAmberDir             } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getCobaltDir            } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { selectCurrentOrExisting } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow COBALT_NORMALISATION {
    take:
    // Sample data
    ch_amber                // channel: [mandatory] [ meta, amber_dir ]
    ch_cobalt               // channel: [mandatory] [ meta, cobalt_dir ]

    // Reference data
    genome_version          // channel: [mandatory] genome version
    gc_profile              // channel: [mandatory] /path/to/gc_profile
    copy_number_percentiles // channel: [mandatory] /path/to/copy_number_percentiles
    target_regions_bed      // channel: [mandatory] /path/to/target_regions_bed

    main:
    // Create process input channel
    // channel: [ [amber_dir, ...], [cobalt_dir, ...] ]
    ch_cobalt_inputs = groupByMeta([
        ch_amber,
        ch_cobalt,
    ])
        .map { meta, amber_dir, cobalt_dir ->
            return [
                selectCurrentOrExisting(amber_dir, getAmberDir(meta)),
                selectCurrentOrExisting(cobalt_dir, getCobaltDir(meta)),
            ]
        }
        .collect(flat: false)
        .map { d -> d.transpose() }

    // Run process
    COBALT_PANEL_NORMALISATION(
        ch_cobalt_inputs,
        genome_version,
        gc_profile,
        copy_number_percentiles,
        target_regions_bed,
    )
}
