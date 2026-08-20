//
// PAVE PON creation prepares the panel-specific small variant artefact resource
//

include { PAVE_PON_PANEL_CREATION } from '../../../modules/local/pave/pon_creation/main'

include { getTumorDnaSageDir      } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getTumorDnaSampleName   } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { selectCurrentOrExisting } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow PAVE_PON_CREATION {
    take:
    // Sample data
    ch_sage_dir_somatic // channel: [mandatory] [ meta, sage_dir ]

    // Reference data
    genome_version      // channel: [mandatory] genome version

    main:
    // Create process input channel
    // channel: [ [sage_vcf, ...], [sage_tbi, ...] ]
    ch_pave_inputs = ch_sage_dir_somatic
        .map { meta, sage_dir ->

            def sage_dir_selected = selectCurrentOrExisting(sage_dir, getTumorDnaSageDir(meta))
            def sage_vcf = sage_dir_selected ? sage_dir_selected.resolve("${getTumorDnaSampleName(meta)}.sage.somatic.vcf.gz") : []
            def sage_tbi = sage_dir_selected ? sage_dir_selected.resolve("${getTumorDnaSampleName(meta)}.sage.somatic.vcf.gz.tbi") : []

            return [sage_vcf, sage_tbi]
        }
        .collect(flat: false)
        .map { d -> d.transpose() }

    // Run process
    PAVE_PON_PANEL_CREATION(
        ch_pave_inputs,
        genome_version,
    )
}
