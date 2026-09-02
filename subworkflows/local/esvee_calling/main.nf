//
// ESVEE detects structural variants, and reports breakends and breakpoints.
//

nextflow.enable.types = true

include { ESVEE } from '../../../modules/local/esvee/main'

include { getNormalReduxDirAlignment  } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { getTumorReduxDirAlignment   } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { getInput                    } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalDnaSample          } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalDnaSampleName      } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSample           } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSampleName       } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasInput                    } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { FileType                    } from '../utils_nfcore_oncoanalyser_pipeline/types_enums'
include { groupByMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { joinMeta                    } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { restoreMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { selectCurrentOrExisting     } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow ESVEE_CALLING {
    take:

    // Sample data
    ch_inputs               : Channel<Map>              // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor      : Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal     : Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, redux_dir ]

    // Reference data
    genome_fasta            : Channel<Path>             // channel: [mandatory] /path/to/genome_fasta
    genome_version          : Channel<String>           // channel: [mandatory] genome version
    genome_fai              : Channel<Path>             // channel: [mandatory] /path/to/genome_fai
    genome_dict             : Channel<Path>             // channel: [mandatory] /path/to/genome_dict
    genome_img              : Channel<Path>             // channel: [mandatory]  /path/to/genome_img
    known_fusions           : Channel<Path>             // channel: [mandatory] /path/to/known_fusions
    pon_breakends           : Channel<Path>             // channel: [mandatory] /path/to/pon_sgl
    pon_breakpoints         : Channel<Path>             // channel: [mandatory] /path/to/pon_sv
    decoy_sequences_image   : Channel<Path>             // channel: [mandatory] /path/to/decoy_sequences_image
    repeatmasker_annotations: Channel<Path>             // channel: [mandatory] /path/to/repeatmasker_annotations
    unmap_regions           : Channel<Path>             // channel: [mandatory] /path/to/unmap_regions
    target_regions_bed      : Channel<Path>?            // channel: [optional]  /path/to/target_regions_bed

    // Params
    sequencing_platform     : String                    // string:  [mandatory] sequencing platform

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
            ch_inputs_sorted.runnable_to.map { d -> d + [null, null] },
        )
        .map { meta, tumor_aln, tumor_idx, normal_aln, normal_idx ->

            def meta_esvee = record(
                key: meta.case_id,
                id: meta.case_id,
                tumor_id: getTumorDnaSampleName(meta),
                normal_id: normal_aln ? getNormalDnaSampleName(meta) : null,
            )

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
            ch_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    emit:
    esvee_dir = ch_outputs // channel: [ meta, esvee_dir ]
}
