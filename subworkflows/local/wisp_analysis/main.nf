//
// WISP estimates tumor purity in longitudinal samples using WGS data of the primary
//

nextflow.enable.types = true

include { WISP  } from '../../../modules/local/wisp/main'

include { FileType                    } from '../utils_nfcore_oncoanalyser_pipeline/types'
include { groupByMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { joinMeta                    } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { restoreMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { getInput                    } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getNormalDnaSample          } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getNormalReduxDirAlignment  } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSample           } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSampleName       } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { selectCurrentOrExisting     } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow WISP_ANALYSIS {
    take:
    // Sample data
    ch_inputs: Channel<Map>                  // channel: [mandatory] [ meta ]
    ch_redux_dir: Channel<Tuple<Map, Path>>               // channel: [mandatory] [ meta, redux_dir ]
    ch_amber_dir: Channel<Tuple<Map, Path>>               // channel: [mandatory] [ meta, amber_dir ]
    ch_cobalt_dir: Channel<Tuple<Map, Path>>              // channel: [mandatory] [ meta, cobalt_dir ]
    ch_sage_append_dir_somatic: Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, sage_append_dir ]

    // Reference data
    genome_fasta: Channel<Path>               // channel: [mandatory] /path/to/genome_fasta
    genome_fai: Channel<Path>                 // channel: [mandatory] /path/to/genome_fai

    // Params
    targeted_mode: Boolean              // boolean: [mandatory] Set targeted mode

    main:
    // Select input sources then sort
    // channel: runnable: [ meta, purple_dir (primary), amber_dir (primary), normal_aln (primary), redux_dir (longitudinal), amber_dir (longitudinal), cobalt_dir (longitudinal), sage_append_dir (longitudinal) ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = groupByMeta([
        ch_redux_dir,
        ch_amber_dir,
        ch_cobalt_dir,
        ch_sage_append_dir_somatic,
    ])
        .map { meta, longitudinal_redux_dir, longitudinal_amber_dir, longitudinal_cobalt_dir, longitudinal_sage_append_dir ->

            def primary_normal_redux_dir = getInput(getNormalDnaSample(meta), FileType.REDUX_DIR)
            def (primary_normal_aln, _primary_normal_idx) = getNormalReduxDirAlignment(meta, primary_normal_redux_dir)

            def primary_purple_dir = getInput(getTumorDnaSample(meta), FileType.PURPLE_DIR)
            def primary_amber_dir = getInput(getTumorDnaSample(meta), FileType.AMBER_DIR)

            def longitudinal_redux_dir_selected = selectCurrentOrExisting(longitudinal_redux_dir, getInput(getTumorDnaSample(meta), FileType.REDUX_DIR))

            return [
              meta,
              primary_purple_dir,
              primary_amber_dir,
              primary_normal_aln,
              longitudinal_redux_dir_selected,
              longitudinal_amber_dir,
              longitudinal_cobalt_dir,
              longitudinal_sage_append_dir,
            ]
        }
        .branch { meta, primary_purple_dir, primary_amber_dir, primary_normal_aln, longitudinal_redux_dir, longitudinal_amber_dir, longitudinal_cobalt_dir, longitudinal_sage_append_dir ->

            def runnable
            if (targeted_mode) {
                runnable = primary_purple_dir && longitudinal_sage_append_dir
            } else {
                runnable = primary_purple_dir && longitudinal_sage_append_dir && longitudinal_cobalt_dir
            }

            runnable: runnable
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_wisp, purple_dir (primary), amber_dir (primary), normal_aln (primary), redux_dir (longitudinal), amber_dir (longitudinal), cobalt_dir (longitudinal), sage_append_dir (longitudinal) ]
    ch_wisp_inputs = ch_inputs_sorted.runnable
        .map { d ->

            def meta = d[0]
            def inputs = d[1..-1]

            def meta_wisp = [
                key: meta.case_id,
                id: meta.case_id,
                patient_id: meta.patient_id,
                primary_id: getTumorDnaSampleName(meta, primary: true),
                longitudinal_id: getTumorDnaSampleName(meta, primary: false),
            ]

            return [meta_wisp] + inputs
        }

    // Run process
    WISP(
        ch_wisp_inputs,
        genome_fasta,
        genome_fai,
        targeted_mode,
    )

    // channel: [ meta, wisp_dir ]
    emit:
    wisp_dir = restoreMeta(channel.topic('wisp_dir'), ch_inputs)
}
