//
// PEACH infers germline haplotypes and reports relevant pharmacogenomics
//

nextflow.enable.types = true

include { PEACH  } from '../../../modules/local/peach/main'

include { FileType                 } from '../utils_nfcore_oncoanalyser_pipeline/types_enums'
include { groupByMeta              } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { joinMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { restoreMeta              } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { getInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalDnaSample       } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalDnaSampleName   } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSample        } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSampleName    } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getPurpleGermlineVcf   } from '../utils_nfcore_oncoanalyser_pipeline/accessors_outputs'
include { hasInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasNormalDna             } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { selectCurrentOrExisting  } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow PEACH_CALLING {
    take:
    // Sample data
    ch_inputs: Channel<Map>                 // channel: [mandatory] [ meta ]
    ch_purple_dir: Channel<Tuple<Map, Path>>             // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    peach_haplotypes: Channel<Path>          // channel: [mandatory] /path/to/peach_haplotypes
    peach_haplotype_functions: Channel<Path> // channel: [mandatory] /path/to/peach_haplotype_functions
    peach_drug_info: Channel<Path>           // channel: [mandatory] /path/to/peach_drug_info

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
            def has_smlv_vcf = getPurpleGermlineVcf(tumor_id, purple_dir)?.exists() ?: false

            runnable: has_smlv_vcf && has_normal && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_peach, purple_germline_smlv_vcf ]
    ch_peach_inputs = ch_inputs_sorted.runnable
        .map { meta, purple_dir ->

            def meta_peach = record(
                key: meta.case_id,
                id: meta.case_id,
                sample_id: getNormalDnaSampleName(meta),
            )

            def purple_germline_smlv_vcf = getPurpleGermlineVcf(getTumorDnaSampleName(meta), purple_dir)[0]

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
