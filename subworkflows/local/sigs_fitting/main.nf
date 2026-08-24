//
// Sigs fits trinucleotide signature definitions with sample SNV counts
//

nextflow.enable.types = true

include { SIGS  } from '../../../modules/local/sigs/main'

include { FileType                 } from '../utils_nfcore_oncoanalyser_pipeline/types_enums'
include { groupByMeta              } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { joinMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { restoreMeta              } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { getInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSample        } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSampleName    } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getPurpleSomaticVcf    } from '../utils_nfcore_oncoanalyser_pipeline/accessors_outputs'
include { hasInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasNormalDna             } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasTumorDna              } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { selectCurrentOrExisting  } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow SIGS_FITTING {
    take:
    // Sample data
    ch_inputs: Channel<Map>       // channel: [mandatory] [ meta ]
    ch_purple_dir: Channel<Tuple<Map, Path>>   // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    sigs_signatures: Channel<Path> // channel: [mandatory] /path/to/sigs_signatures

    main:
    // Select input sources then sort
    // channel: runnable: [ meta, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_purple_dir
        .map { meta, purple_dir ->
            return [meta, selectCurrentOrExisting(purple_dir, getInput(getTumorDnaSample(meta), FileType.PURPLE_DIR))]
        }
        .branch { meta, purple_dir ->

            def has_tumor_normal_dna = hasTumorDna(meta) && hasNormalDna(meta)

            def has_smlv_vcf = false
            if (has_tumor_normal_dna && purple_dir) {
                def tumor_id = getTumorDnaSampleName(meta)
                has_smlv_vcf = getPurpleSomaticVcf(tumor_id, purple_dir).exists()
            }

            def has_existing = hasInput(getTumorDnaSample(meta), FileType.SIGS_DIR)

            runnable: has_tumor_normal_dna && purple_dir && has_smlv_vcf && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_sigs, smlv_vcf ]
    ch_sigs_inputs = ch_inputs_sorted.runnable
        .map { meta, purple_dir ->

            def tumor_id = getTumorDnaSampleName(meta)

            def meta_sigs = record(
                key: meta.case_id,
                id: meta.case_id,
                sample_id: tumor_id,
            )

            def smlv_vcf = getPurpleSomaticVcf(tumor_id, purple_dir)

            return [meta_sigs, smlv_vcf]
        }

    // Run process
    SIGS(
        ch_sigs_inputs,
        sigs_signatures,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, sigs_dir ]
    ch_outputs = channel.empty()
        .mix(
            restoreMeta(channel.topic('sigs_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    emit:
    sigs_dir = ch_outputs // channel: [ meta, sigs_dir ]
}
