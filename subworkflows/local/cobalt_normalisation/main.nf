//
// COBALT normalisation prepares the panel-specific target region normalisation resource
//

include { COBALT_PANEL_NORMALISATION } from '../../../modules/local/cobalt/panel_normalisation/main'

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
        copy_number_percentiles,
        target_regions_bed,
    )
}
