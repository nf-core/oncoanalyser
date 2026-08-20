//
// Qsee calculates and visualises QC metrics
//

include { QSEE } from '../../../modules/local/qsee/main'
include { groupByMeta; joinMeta; restoreMeta } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { getCobaltDir; getEsveeDir; getNormalDnaBamtoolsDir; getNormalDnaReduxDir; getNormalDnaSampleName; getNormalReduxTsvs; getPurpleDir; getTumorDnaBamtoolsDir; getTumorDnaReduxDir; getTumorDnaSampleName; getTumorReduxTsvs; selectCurrentOrExisting } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow QSEE_METRICS {
    take:
    // Sample data
    ch_inputs                // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor       // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal      // channel: [mandatory] [ meta, redux_dir ]
    ch_bamtools_dir_tumor    // channel: [mandatory] [ meta, bamtools_dir ]
    ch_bamtools_dir_normal   // channel: [optional]  [ meta, bamtools_dir ]
    ch_cobalt_dir            // channel: [optional]  [ meta, cobalt_dir ]
    ch_esvee_dir             // channel: [optional]  [ meta, esvee_dir ]
    ch_purple_dir            // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    driver_gene_panel        // channel: [mandatory] /path/to/driver_gene_panel
    qsee_cohort_percentiles  // channel: [mandatory] /path/to/cohort_percentiles

    // Params
    sequencing_platform      // string:  [mandatory] sequencing platform
    targeted_mode            // boolean: [mandatory] Set targeted mode

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

            def redux_tumor_dir_selected = selectCurrentOrExisting(redux_dir_tumor, getTumorDnaReduxDir(meta))
            def redux_normal_dir_selected = selectCurrentOrExisting(redux_dir_normal, getNormalDnaReduxDir(meta))

            def redux_tsvs_tumor = getTumorReduxTsvs(meta, redux_tumor_dir_selected)
            def redux_tsvs_normal = getNormalReduxTsvs(meta, redux_normal_dir_selected)

            return [
                meta,
                redux_tsvs_tumor,
                redux_tsvs_normal,
                selectCurrentOrExisting(bamtools_dir_tumor, getTumorDnaBamtoolsDir(meta)),
                selectCurrentOrExisting(bamtools_dir_normal, getNormalDnaBamtoolsDir(meta)),
                selectCurrentOrExisting(cobalt_dir, getCobaltDir(meta)),
                selectCurrentOrExisting(esvee_dir, getEsveeDir(meta)),
                selectCurrentOrExisting(purple_dir, getPurpleDir(meta)),
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

            def meta_qsee = [
                key: meta.case_id,
                id: meta.case_id,
                tumor_id: getTumorDnaSampleName(meta),
            ]

            if (redux_tsvs_normal || bamtools_dir_normal) {
                meta_qsee.normal_id = getNormalDnaSampleName(meta)
            }

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
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    qsee_dir = ch_outputs // channel: [ meta, qsee_dir ]
}
