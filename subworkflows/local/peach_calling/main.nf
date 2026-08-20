//
// PEACH infers germline haplotypes and reports relevant pharmacogenomics
//

include { PEACH } from '../../../modules/local/peach/main'
include { groupByMeta             } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { joinMeta                } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { restoreMeta             } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { getNormalDnaSampleName  } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getPurpleDir            } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getTumorDnaSampleName   } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { hasNormalDna            } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { hasPeachDir             } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { selectCurrentOrExisting } from '../utils_nfcore_oncoanalyser_pipeline/utils'

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
                selectCurrentOrExisting(purple_dir, getPurpleDir(meta)),
            ]
        }
        .branch { meta, purple_dir ->

            def has_normal = hasNormalDna(meta)
            def has_existing = hasPeachDir(meta)

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
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    peach_dir = ch_outputs // channel: [ meta, peach_dir ]
}
