//
// COBALT calculates read ratios between tumor and normal samples
//

nextflow.enable.types = true

include { COBALT } from '../../../modules/local/cobalt/run/main'

include { getNormalReduxDirAlignment } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { getTumorReduxDirAlignment  } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { getInput                   } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getLongitudinalSampleName  } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalDnaSample         } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalDnaSampleName     } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSample          } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSampleName      } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasInput                   } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { groupByMeta                } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { joinMeta                   } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { restoreMeta                } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { FileType                   } from '../utils_nfcore_oncoanalyser_pipeline/types_enums'
include { selectCurrentOrExisting    } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow COBALT_PROFILING {
    take:
    // Sample data
    ch_inputs                   : Channel<Map>              // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor          : Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal         : Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, redux_dir ]

    // Reference data
    genome_fasta                : Channel<Path>             // channel: [mandatory] /path/to/genome_fasta
    genome_version              : Channel<String>           // channel: [mandatory] genome version
    genome_fai                  : Channel<Path>             // channel: [mandatory] /path/to/genome_fai
    gc_profile                  : Channel<Path>             // channel: [mandatory] /path/to/gc_profile
    diploid_bed                 : Channel<Path>?            // channel: [optional]  /path/to/diploid_bed
    target_regions_normalisation: Channel<Path>?            // channel: [optional]  /path/to/target_regions_normalisation
    targeted_mode               : Boolean                   // boolean: [mandatory] Set targeted mode
    purity_estimate_mode        : Boolean                   // boolean: [mandatory] Set purity estimate mode

    main:
    // Select input sources then sort
    // NOTE(SW): germline mode is not currently supported
    // channel: runnable: [ meta, tumor_aln, tumor_idx, normal_aln, normal_idx ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = groupByMeta([
        ch_redux_dir_tumor,
        ch_redux_dir_normal,
    ])
        .map { meta, redux_dir_tumor, redux_dir_normal ->

            def redux_dir_tumor_selected = selectCurrentOrExisting(redux_dir_tumor, getInput(getTumorDnaSample(meta), FileType.REDUX_DIR))
            def redux_dir_normal_selected = selectCurrentOrExisting(redux_dir_normal, getInput(getNormalDnaSample(meta), FileType.REDUX_DIR))

            def (tumor_aln, tumor_idx) = getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)
            def (normal_aln, normal_idx) = getNormalReduxDirAlignment(meta, redux_dir_normal_selected)

            return [meta, tumor_aln, tumor_idx, normal_aln, normal_idx]

        }
        .branch { meta, tumor_aln, tumor_idx, normal_aln, normal_idx ->
            def has_existing = hasInput(getTumorDnaSample(meta), FileType.COBALT_DIR)
            runnable_tn: tumor_aln && normal_aln && ! has_existing
            runnable_to: tumor_aln && ! has_existing
            skip: true
                return meta
        }

    // First set diploid BED input for tumor/normal and tumor only samples
    // NOTE(SW): since the diploid BED is provided as a channel, I seem to be only able to include via channel ops
    // channel: [ meta, tumor_aln, tumor_idx, normal_aln, normal_idx, diploid_bed ]
    ch_inputs_runnable = channel.empty()
        .mix(
            ch_inputs_sorted.runnable_tn.map { d -> d + [null] },
            ch_inputs_sorted.runnable_to.combine(diploid_bed),
        )

    // Create process input channel
    // channel: sample_data: [ meta_cobalt, tumor_aln, tumor_idx, normal_aln, normal_idx ]
    // channel: diploid_bed: [ diploid_bed ]
    ch_cobalt_inputs = ch_inputs_runnable
        .multiMap { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, ref_diploid_bed ->

            def tumor_id
            if (purity_estimate_mode) {
                tumor_id = getLongitudinalSampleName(meta)
            } else {
                tumor_id = getTumorDnaSampleName(meta)
            }

            def meta_cobalt = record(
                key: meta.case_id,
                id: meta.case_id,
                tumor_id: tumor_id,
                normal_id: normal_aln ? getNormalDnaSampleName(meta) : null,
            )

            sample_data: [meta_cobalt, tumor_aln, tumor_idx, normal_aln, normal_idx]
            diploid_bed: ref_diploid_bed
        }

    // Run process
    COBALT(
        ch_cobalt_inputs.sample_data,
        genome_fasta,
        genome_version,
        genome_fai,
        gc_profile,
        ch_cobalt_inputs.diploid_bed,
        target_regions_normalisation,
        targeted_mode,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, cobalt_dir ]
    ch_outputs = channel.empty()
        .mix(
            restoreMeta(channel.topic('cobalt_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    emit:
    cobalt_dir = ch_outputs // channel: [ meta, cobalt_dir ]
}
