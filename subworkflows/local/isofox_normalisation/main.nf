//
// ISOFOX normalisation prepares panel-specific TPM normalisation resource
//

include { ISOFOX_PANEL_NORMALISATION  } from '../../../modules/local/isofox/panel_normalisation/main'

include { FileType                 } from '../utils_nfcore_oncoanalyser_pipeline/types'
include { getInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorRnaSample        } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { selectCurrentOrExisting  } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow ISOFOX_NORMALISATION {
    take:
    // Sample data
    ch_isofox                // channel: [mandatory] [ meta, isofox_dir ]

    // Reference data
    genome_version           // channel: [mandatory] genome version
    isofox_gene_ids          // channel: [mandatory]  /path/to/gene_ids
    isofox_gene_distribution // channel: [mandatory] /path/to/isofox_gene_distribution

    main:
    // Create process input channel
    // channel: [ [isofox_dir, ...] ]
    ch_isofox_inputs = ch_isofox
        .map { meta, isofox_dir ->
            return selectCurrentOrExisting(isofox_dir, getInput(getTumorRnaSample(meta), FileType.ISOFOX_DIR))
        }
        .collect()

    // Run process
    ISOFOX_PANEL_NORMALISATION(
        ch_isofox_inputs,
        genome_version,
        isofox_gene_ids,
        isofox_gene_distribution,
    )

    emit:
    isofox_normalisation_csv = channel.topic('isofox_normalisation_csv')
}
