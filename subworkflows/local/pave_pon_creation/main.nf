//
// PAVE PON creation prepares the panel-specific small variant artefact resource
//

nextflow.enable.types = true

include { PAVE_PON_PANEL_CREATION } from '../../../modules/local/pave/pon_creation/main'

include { getInput                } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSample       } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSampleName   } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getSageSomaticVcf       } from '../utils_nfcore_oncoanalyser_pipeline/accessors_outputs'
include { getVcfTbi               } from '../utils_nfcore_oncoanalyser_pipeline/accessors_outputs'
include { FileType                } from '../utils_nfcore_oncoanalyser_pipeline/types_enums'
include { selectCurrentOrExisting } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow PAVE_PON_CREATION {
    take:
    // Sample data
    ch_sage_dir_somatic: Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, sage_dir ]

    // Reference data
    genome_version     : Channel<String>           // channel: [mandatory] genome version

    main:
    // Create process input channel
    // channel: [ [sage_vcf, ...], [sage_tbi, ...] ]
    ch_pave_inputs = ch_sage_dir_somatic
        .map { meta, sage_dir ->

            def sage_dir_selected = selectCurrentOrExisting(sage_dir, getInput(getTumorDnaSample(meta), FileType.SAGE_DIR))
            def sage_vcf = getSageSomaticVcf(getTumorDnaSampleName(meta), sage_dir_selected)
            def sage_tbi = getVcfTbi(sage_vcf)

            return [sage_vcf, sage_tbi]
        }
        .collect(flat: false)
        .map { d -> d.transpose() }

    // Run process
    PAVE_PON_PANEL_CREATION(
        ch_pave_inputs,
        genome_version,
    )

    // Set outputs
    // channel: [ pave_pon_artefacts ]
    ch_outputs = channel.topic('pave_pon_panel_creation_artefacts')

    emit:
    pave_pon_artefacts = ch_outputs // channel: [ pave_pon_artefacts ]
}
