//
// TEAL performs characterisation of telomeric features and rearrangements
//

include { TEAL_PREP  } from '../../../modules/local/teal/prep/main'
include { TEAL_PIPELINE  } from '../../../modules/local/teal/pipeline/main'

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
include { selectCurrentOrExisting     } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow TEAL_CHARACTERISATION {
    take:
    // Sample data
    ch_inputs              // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor     // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal    // channel: [mandatory] [ meta, redux_dir ]
    ch_bamtools_dir_tumor  // channel: [mandatory] [ meta, bamtools_dir ]
    ch_bamtools_dir_normal // channel: [mandatory] [ meta, bamtools_dir ]
    ch_cobalt_dir          // channel: [mandatory] [ meta, cobalt_dir ]
    ch_purple_dir          // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_fasta           // channel: [mandatory] /path/to/genome_fasta
    genome_version         // channel: [mandatory] genome version
    genome_fai             // channel: [mandatory] /path/to/genome_fai

    // Params
    sequencing_platform    // string:  [mandatory] sequencing platform

    main:
    //
    // STEP: Handle inputs
    //
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
            runnable: tumor_aln || normal_aln
            skip: true
                return meta
        }

    //
    // MODULE: TEAL prep
    //
    // Create process input channel
    // channel: [ meta_teal, tumor_aln, tumor_idx, normal_aln, normal_idx ]
    ch_teal_prep_inputs = ch_inputs_sorted.runnable
        .map { meta, tumor_aln, tumor_idx, normal_aln, normal_idx ->

            def meta_teal = [
                key: meta.case_id,
                id: meta.case_id,
            ]

            if (tumor_aln) {
                meta_teal.tumor_id = getTumorDnaSampleName(meta)
            }

            if (normal_aln) {
                meta_teal.normal_id = getNormalDnaSampleName(meta)
            }

            return [meta_teal, tumor_aln, tumor_idx, normal_aln, normal_idx]
        }

    // Run process
    TEAL_PREP(
        ch_teal_prep_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        sequencing_platform,
    )

    // Restore metas
    // channel: [ meta, teal_bam, teal_bai ]
    ch_tumor_teal_bam = restoreMeta(channel.topic('teal_prep_tumor_bam'), ch_inputs)
    ch_normal_teal_bam = restoreMeta(channel.topic('teal_prep_normal_bam'), ch_inputs)

    //
    // MODULE: TEAL pipeline
    //
    // Select input sources then sort
    // channel: runnable: [ meta, teal_bam_tumor, teal_bai_tumor, teal_bam_normal, teal_bai_normal, bamtools_dir_tumor, bamtools_dir_normal, cobalt_dir, purple_dir ]
    // channel: skip: [ meta ]
    ch_teal_pipeline_inputs_sorted = groupByMeta([
        ch_tumor_teal_bam,
        ch_normal_teal_bam,
        ch_bamtools_dir_tumor,
        ch_bamtools_dir_normal,
        ch_cobalt_dir,
        ch_purple_dir,
    ])
        .map { meta, teal_bam_tumor, teal_bai_tumor, teal_bam_normal, teal_bai_normal, bamtools_dir_tumor, bamtools_dir_normal, cobalt_dir, purple_dir ->
            return [
                meta,
                teal_bam_tumor,
                teal_bai_tumor,
                teal_bam_normal,
                teal_bai_normal,
                selectCurrentOrExisting(bamtools_dir_tumor, getInput(getTumorDnaSample(meta), FileType.BAMTOOLS_DIR)),
                selectCurrentOrExisting(bamtools_dir_normal, getInput(getNormalDnaSample(meta), FileType.BAMTOOLS_DIR)),
                selectCurrentOrExisting(cobalt_dir, getInput(getTumorDnaSample(meta), FileType.COBALT_DIR)),
                selectCurrentOrExisting(purple_dir, getInput(getTumorDnaSample(meta), FileType.PURPLE_DIR)),
            ]
        }
        .branch { meta, teal_bam_tumor, teal_bai_tumor, teal_bam_normal, teal_bai_normal, bamtools_dir_tumor, bamtools_dir_normal, cobalt_dir, purple_dir ->

            def has_tumor = teal_bam_tumor && bamtools_dir_tumor && purple_dir
            def has_normal = teal_bam_normal && bamtools_dir_normal

            runnable: (has_tumor || has_normal) && cobalt_dir
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_teal, teal_bam_tumor, teal_bai_tumor, teal_bam_normal, teal_bai_normal, bamtools_dir_tumor, bamtools_dir_normal, cobalt_dir, purple_dir ]
    ch_teal_pipeline_inputs = ch_teal_pipeline_inputs_sorted.runnable
        .map { meta, teal_bam_tumor, teal_bai_tumor, teal_bam_normal, teal_bai_normal, bamtools_dir_tumor, bamtools_dir_normal, cobalt_dir, purple_dir ->

            def meta_teal = [
                key: meta.case_id,
                id: meta.case_id,
            ]

            if (teal_bam_tumor) {
                meta_teal.tumor_id = getTumorDnaSampleName(meta)
            }

            if (teal_bam_normal) {
                meta_teal.normal_id = getNormalDnaSampleName(meta)
            }

            return [meta_teal, teal_bam_tumor, teal_bai_tumor, teal_bam_normal, teal_bai_normal, bamtools_dir_tumor, bamtools_dir_normal, cobalt_dir, purple_dir]
        }

    // Run process
    TEAL_PIPELINE(
        ch_teal_pipeline_inputs,
        genome_version,
        sequencing_platform,
    )

    // channel: [ meta, teal_bam, teal_bai ]
    // channel: [ meta, teal_tsvs ]
    emit:
    tumor_bam = ch_tumor_teal_bam
    normal_bam = ch_normal_teal_bam
    teal_tsvs = restoreMeta(channel.topic('teal_tsvs'), ch_inputs)
}
