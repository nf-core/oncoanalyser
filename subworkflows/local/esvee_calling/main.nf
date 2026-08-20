//
// ESVEE detects structural variants, and reports breakends and breakpoints.
//

include { ESVEE  } from '../../../modules/local/esvee/main'

include { FileType                    } from '../utils_nfcore_oncoanalyser_pipeline/types'
include { groupByMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { joinMeta                    } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { restoreMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { getInput                    } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getNormalDnaSample          } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getNormalDnaSampleName      } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getNormalReduxDirAlignment  } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSample           } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSampleName       } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorReduxDirAlignment   } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { hasInput                    } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { selectCurrentOrExisting     } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow ESVEE_CALLING {
    take:

    // Sample data
    ch_inputs                // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor       // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal      // channel: [mandatory] [ meta, redux_dir ]

    // Reference data
    genome_fasta             // channel: [mandatory] /path/to/genome_fasta
    genome_version           // channel: [mandatory] genome version
    genome_fai               // channel: [mandatory] /path/to/genome_fai
    genome_dict              // channel: [mandatory] /path/to/genome_dict
    genome_img               // channel: [optional]  /path/to/genome_img
    known_fusions            // channel: [mandatory] /path/to/known_fusions
    pon_breakends            // channel: [mandatory] /path/to/pon_sgl
    pon_breakpoints          // channel: [mandatory] /path/to/pon_sv
    decoy_sequences_image    // channel: [mandatory] /path/to/decoy_sequences_image
    repeatmasker_annotations // channel: [mandatory] /path/to/repeatmasker_annotations
    unmap_regions            // channel: [mandatory] /path/to/unmap_regions
    target_regions_bed       // channel: [optional]  /path/to/target_regions_bed

    // Params
    sequencing_platform      // string:  [mandatory] sequencing platform

    main:
    // Select input sources then sort
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
            def has_existing = hasInput(getTumorDnaSample(meta), FileType.ESVEE_DIR)

            runnable_tn: tumor_aln && normal_aln && ! has_existing
            runnable_to: tumor_aln && ! has_existing
                return [meta, tumor_aln, tumor_idx]
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_esvee, tumor_aln, tumor_idx, normal_aln, normal_idx ]
    ch_esvee_inputs = channel.empty()
        .mix(
            ch_inputs_sorted.runnable_tn,
            ch_inputs_sorted.runnable_to.map { d -> d + [[], []] },
        )
        .map { meta, tumor_aln, tumor_idx, normal_aln, normal_idx ->

            def meta_esvee = [
                key: meta.case_id,
                id: meta.case_id,
                tumor_id: getTumorDnaSampleName(meta),
            ]

            if (normal_aln) {
                meta_esvee.normal_id = getNormalDnaSampleName(meta)
            }

            return [meta_esvee, tumor_aln, tumor_idx, normal_aln, normal_idx]
        }

    // Run process
    ESVEE(
        ch_esvee_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        genome_img,
        pon_breakends,
        pon_breakpoints,
        decoy_sequences_image,
        known_fusions,
        repeatmasker_annotations,
        unmap_regions,
        target_regions_bed,
        sequencing_platform,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, esvee_dir ]
    ch_outputs = channel.empty()
        .mix(
            restoreMeta(channel.topic('esvee_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    esvee_dir = ch_outputs // channel: [ meta, esvee_dir ]
}
