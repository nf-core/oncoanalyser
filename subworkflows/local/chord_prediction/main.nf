//
// CHORD predicts HR status for tumor samples
//

nextflow.enable.types = true

include { CHORD  } from '../../../modules/local/chord/main'

include { FileType                 } from '../utils_nfcore_oncoanalyser_pipeline/types'
include { groupByMeta              } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { joinMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { restoreMeta              } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { getInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSample        } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSampleName    } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { hasInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { hasNormalDna             } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { hasTumorDna              } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { selectCurrentOrExisting  } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow CHORD_PREDICTION {
    take:
    // Sample data
    ch_inputs: Channel<Map>      // channel: [mandatory] [ meta ]
    ch_purple_dir: Channel<Tuple<Map, Path>>  // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_fasta: Channel<Path>   // channel: [mandatory] /path/to/genome_fasta
    genome_fai: Channel<Path>     // channel: [mandatory] /path/to/genome_fai
    genome_dict: Channel<Path>    // channel: [mandatory] /path/to/genome_dict

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
            def has_sv_vcf = false
            if (has_tumor_normal_dna && purple_dir) {
                def tumor_id = getTumorDnaSampleName(meta)
                has_smlv_vcf = purple_dir.resolve("${tumor_id}.purple.somatic.vcf.gz").exists()
                has_sv_vcf = purple_dir.resolve("${tumor_id}.purple.sv.vcf.gz").exists()
            }

            def has_existing = hasInput(getTumorDnaSample(meta), FileType.CHORD_DIR)

            runnable: has_tumor_normal_dna && purple_dir && has_smlv_vcf && has_sv_vcf && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_chord, smlv_vcf, sv_vcf ]
    ch_chord_inputs = ch_inputs_sorted.runnable
        .map { meta, purple_dir ->

            def tumor_id = getTumorDnaSampleName(meta)

            def meta_chord = record(
                key: meta.case_id,
                id: meta.case_id,
                sample_id: tumor_id,
            )

            def smlv_vcf = purple_dir.resolve("${tumor_id}.purple.somatic.vcf.gz")
            def sv_vcf = purple_dir.resolve("${tumor_id}.purple.sv.vcf.gz")

            return [meta_chord, smlv_vcf, sv_vcf]
        }

    // Run process
    CHORD(
        ch_chord_inputs,
        genome_fasta,
        genome_fai,
        genome_dict,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, chord_dir ]
    ch_outputs = channel.empty()
        .mix(
            restoreMeta(channel.topic('chord_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    emit:
    chord_dir = ch_outputs // channel: [ meta, chord_dir ]
}
