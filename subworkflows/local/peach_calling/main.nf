//
// PEACH infers germline haplotypes and reports relevant pharmacogenomics
//

include { PEACH  } from '../../../modules/local/peach/main'

include { FileType                 } from '../utils_nfcore_oncoanalyser_pipeline/types'
include { groupByMeta              } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { joinMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { restoreMeta              } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { getInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getNormalDnaSample       } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getNormalDnaSampleName   } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSample        } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSampleName    } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { hasInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { hasNormalDna             } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { selectCurrentOrExisting  } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow PEACH_CALLING {
    take:
    // Sample data
    ch_inputs                 // channel: [mandatory] [ meta ]
    ch_purple_dir             // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    peach_haplotypes          // channel: [mandatory] /path/to/peach_haplotypes
    peach_haplotype_functions // channel: [mandatory] /path/to/peach_haplotype_functions
    peach_drug_info           // channel: [mandatory] /path/to/peach_drug_info

    main:
    // Select input sources then sort
    // channel: runnable: [ meta, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_purple_dir
        .map { meta, purple_dir ->
            return [
                meta,
                selectCurrentOrExisting(purple_dir, getInput(getTumorDnaSample(meta), FileType.PURPLE_DIR)),
            ]
        }
        .branch { meta, purple_dir ->

            def has_normal = hasNormalDna(meta)
            def has_existing = hasInput(getNormalDnaSample(meta), FileType.PEACH_DIR)

            def tumor_id = getTumorDnaSampleName(meta)
            def has_smlv_vcf = purple_dir ? purple_dir.resolve("${tumor_id}.purple.germline.vcf.gz").exists() : false

            runnable: has_smlv_vcf && has_normal && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_peach, purple_germline_smlv_vcf ]
    ch_peach_inputs = ch_inputs_sorted.runnable
        .map { meta, purple_dir ->

            def meta_peach = [
                key: meta.case_id,
                id: meta.case_id,
                sample_id: getNormalDnaSampleName(meta),
            ]

            def purple_germline_smlv_vcf = purple_dir.resolve("${getTumorDnaSampleName(meta)}.purple.germline.vcf.gz")

            return [meta_peach, purple_germline_smlv_vcf]
        }

    // Run process
    PEACH(
        ch_peach_inputs,
        peach_haplotypes,
        peach_haplotype_functions,
        peach_drug_info,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, peach_dir ]
    ch_outputs = channel.empty()
        .mix(
            restoreMeta(channel.topic('peach_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    emit:
    peach_dir = ch_outputs // channel: [ meta, peach_dir ]
}
