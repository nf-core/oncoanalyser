//
// ISOFOX normalisation prepares panel-specific TPM normalisation resource
//

nextflow.enable.types = true

include { ISOFOX_PANEL_NORMALISATION } from '../../../modules/local/isofox/panel_normalisation/main'

include { getInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorRnaSample        } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { FileType                 } from '../utils_nfcore_oncoanalyser_pipeline/types_enums'
include { selectCurrentOrExisting  } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow ISOFOX_NORMALISATION {
    take:
    // Sample data
    ch_isofox: Channel<Tuple<Map, Path>>                // channel: [mandatory] [ meta, isofox_dir ]

    // Reference data
    genome_version: Channel<String>           // channel: [mandatory] genome version
    isofox_gene_ids: Channel<Path>          // channel: [mandatory]  /path/to/gene_ids
    isofox_gene_distribution: Channel<Path> // channel: [mandatory] /path/to/isofox_gene_distribution

    main:
    // Create process input channel
    // channel: [ [isofox_dir, ...] ]
    ch_isofox_inputs = ch_isofox
        .map { meta, isofox_dir ->
            return selectCurrentOrExisting(isofox_dir, getInput(getTumorRnaSample(meta), FileType.ISOFOX_DIR))
        }
        .filter { it != null }
        .toList()

    // Run process
    ISOFOX_PANEL_NORMALISATION(
        ch_isofox_inputs,
        genome_version,
        isofox_gene_ids,
        isofox_gene_distribution,
    )

    // Set outputs
    // channel: [ meta, isofox_normalisation_csv ]
    ch_outputs = channel.topic('isofox_normalisation_csv')

    emit:
    isofox_normalisation_csv = ch_outputs // channel: [ meta, isofox_normalisation_csv ]
}
