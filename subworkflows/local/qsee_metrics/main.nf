//
// Qsee calculates and visualises QC metrics
//

nextflow.enable.types = true

include { QSEE } from '../../../modules/local/qsee/main'

include { getNormalReduxTsvs      } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { getTumorReduxTsvs       } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { getInput                } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalDnaSample      } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalDnaSampleName  } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSample       } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSampleName   } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { FileType                } from '../utils_nfcore_oncoanalyser_pipeline/types_enums'
include { groupByMeta             } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { joinMeta                } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { restoreMeta             } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { selectCurrentOrExisting } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow QSEE_METRICS {
    take:
    // Sample data
    ch_inputs: Channel<Map>                // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor: Channel<Tuple<Map, Path>>       // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal: Channel<Tuple<Map, Path>>      // channel: [mandatory] [ meta, redux_dir ]
    ch_bamtools_dir_tumor: Channel<Tuple<Map, Path>>    // channel: [mandatory] [ meta, bamtools_dir ]
    ch_bamtools_dir_normal: Channel<Tuple<Map, Path>>   // channel: [optional]  [ meta, bamtools_dir ]
    ch_cobalt_dir: Channel<Tuple<Map, Path>>            // channel: [optional]  [ meta, cobalt_dir ]
    ch_esvee_dir: Channel<Tuple<Map, Path>>             // channel: [optional]  [ meta, esvee_dir ]
    ch_purple_dir: Channel<Tuple<Map, Path>>            // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    driver_gene_panel: Channel<Path>        // channel: [mandatory] /path/to/driver_gene_panel
    qsee_cohort_percentiles: Channel<Path>  // channel: [mandatory] /path/to/cohort_percentiles

    // Params
    sequencing_platform: String      // string:  [mandatory] sequencing platform
    targeted_mode: Boolean            // boolean: [mandatory] Set targeted mode

    main:
    // Select input sources then sort
    // channel: runnable: [ meta, [redux_tsv_tumor, ...], [redux_tsv_normal, ...], bamtools_tumor_dir, bamtools_normal_dir, cobalt_dir, esvee_dir, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = groupByMeta([
        ch_redux_dir_tumor,
        ch_redux_dir_normal,
        ch_bamtools_dir_tumor,
        ch_bamtools_dir_normal,
        ch_cobalt_dir,
        ch_esvee_dir,
        ch_purple_dir,
    ])
        .map { meta, redux_dir_tumor, redux_dir_normal, bamtools_dir_tumor, bamtools_dir_normal, cobalt_dir, esvee_dir, purple_dir ->

            def redux_tumor_dir_selected = selectCurrentOrExisting(redux_dir_tumor, getInput(getTumorDnaSample(meta), FileType.REDUX_DIR))
            def redux_normal_dir_selected = selectCurrentOrExisting(redux_dir_normal, getInput(getNormalDnaSample(meta), FileType.REDUX_DIR))

            def redux_tsvs_tumor = getTumorReduxTsvs(meta, redux_tumor_dir_selected)
            def redux_tsvs_normal = getNormalReduxTsvs(meta, redux_normal_dir_selected)

            return [
                meta,
                redux_tsvs_tumor,
                redux_tsvs_normal,
                selectCurrentOrExisting(bamtools_dir_tumor, getInput(getTumorDnaSample(meta), FileType.BAMTOOLS_DIR)),
                selectCurrentOrExisting(bamtools_dir_normal, getInput(getNormalDnaSample(meta), FileType.BAMTOOLS_DIR)),
                selectCurrentOrExisting(cobalt_dir, getInput(getTumorDnaSample(meta), FileType.COBALT_DIR)),
                selectCurrentOrExisting(esvee_dir, getInput(getTumorDnaSample(meta), FileType.ESVEE_DIR)),
                selectCurrentOrExisting(purple_dir, getInput(getTumorDnaSample(meta), FileType.PURPLE_DIR)),
            ]

        }
        .branch { meta, redux_tsvs_tumor, redux_tsvs_normal, bamtools_dir_tumor, bamtools_dir_normal, cobalt_dir, esvee_dir, purple_dir ->

            runnable: bamtools_dir_tumor && purple_dir
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_qsee, [redux_tsv_tumor, ...], [redux_tsv_normal, ...], bamtools_dir_tumor, bamtools_dir_normal, cobalt_dir, esvee_dir, purple_dir ]
    ch_qsee_inputs = ch_inputs_sorted.runnable
        .map { meta, redux_tsvs_tumor, redux_tsvs_normal, bamtools_dir_tumor, bamtools_dir_normal, cobalt_dir, esvee_dir, purple_dir ->

            def meta_qsee = record(
                key: meta.case_id,
                id: meta.case_id,
                tumor_id: getTumorDnaSampleName(meta),
                normal_id: (redux_tsvs_normal || bamtools_dir_normal) ? getNormalDnaSampleName(meta) : null,
            )

            return [meta_qsee, redux_tsvs_tumor, redux_tsvs_normal, bamtools_dir_tumor, bamtools_dir_normal, cobalt_dir, esvee_dir, purple_dir]

        }

    // Run process
    QSEE(
        ch_qsee_inputs,
        driver_gene_panel,
        qsee_cohort_percentiles,
        sequencing_platform,
        targeted_mode,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, qsee_dir ]
    ch_outputs = channel.empty()
        .mix(
            restoreMeta(channel.topic('qsee_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    emit:
    qsee_dir = ch_outputs // channel: [ meta, qsee_dir ]
}
