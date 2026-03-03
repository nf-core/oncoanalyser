//
// Qsee calculates and visualises QC metrics
//

include { QSEE } from '../../../modules/local/qsee/main'

workflow QSEE_METRICS {
    take:
    // Sample data
    ch_inputs                // channel: [mandatory] [ meta ]
    ch_redux_tumor_tsv       // channel: [mandatory] [ meta, redux_tsv, ... ]
    ch_redux_normal_tsv      // channel: [mandatory] [ meta, redux_tsv, ... ]
    ch_bamtools_tumor        // channel: [mandatory] [ meta, metrics_dir ]
    ch_bamtools_normal       // channel: [optional]  [ meta, metrics_dir ]
    ch_cobalt                // channel: [optional]  [ meta, cobalt_dir ]
    ch_esvee                 // channel: [optional]  [ meta, esvee_dir ]
    ch_purple                // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    driver_gene_panel        // channel: [mandatory] /path/to/driver_gene_panel
    qsee_cohort_percentiles  // channel: [mandatory] /path/to/cohort_percentiles

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select and sort inputs
    // channel: [ meta, redux_tumor_tsv, redux_normal_tsv, bamtools_tumor_dir, bamtools_normal_dir, cobalt_dir, esvee_dir, purple_dir ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_redux_tumor_tsv,
        ch_redux_normal_tsv,
        ch_bamtools_tumor,
        ch_bamtools_normal,
        ch_cobalt,
        ch_esvee,
        ch_purple,
    )
        .map { meta,
            tumor_bqr_tsv , tumor_jitter_tsv , tumor_ms_tsv,
            normal_bqr_tsv, normal_jitter_tsv, normal_ms_tsv,
            bamtools_tumor_dir, bamtools_normal_dir,
            cobalt_dir,
            esvee_dir,
            purple_dir ->

            redux_tumor_tsv = [
                tumor_bqr_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_BQR_TSV_TUMOR),
                tumor_dup_freq_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_DUP_FREQ_TSV_TUMOR),
                tumor_ms_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_MS_TSV_TUMOR),
            ]
            redux_tumor_tsv = redux_tumor_tsv.findAll { it != [] }

            redux_normal_tsv = [
                normal_bqr_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_BQR_TSV_NORMAL),
                normal_dup_freq_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_DUP_FREQ_TSV_NORMAL),
                normal_ms_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_MS_TSV_NORMAL),
            ]
            redux_normal_tsv = redux_normal_tsv.findAll { it != [] }

            bamtools_tumor_dir = Utils.selectCurrentOrExisting(bamtools_tumor_dir, meta, Constants.INPUT.BAMTOOLS_DIR_TUMOR)
            bamtools_normal_dir = Utils.selectCurrentOrExisting(bamtools_normal_dir, meta, Constants.INPUT.BAMTOOLS_DIR_NORMAL)

            cobalt_dir = Utils.selectCurrentOrExisting(cobalt_dir, meta, Constants.INPUT.COBALT_DIR)
            esvee_dir = Utils.selectCurrentOrExisting(esvee_dir, meta, Constants.INPUT.ESVEE_DIR)
            purple_dir = Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR)

            return [ meta, redux_tumor_tsv, redux_normal_tsv, bamtools_tumor_dir, bamtools_normal_dir, cobalt_dir, esvee_dir, purple_dir ]
        }
        .branch { meta, redux_tumor_tsv, redux_normal_tsv, bamtools_tumor_dir, bamtools_normal_dir, cobalt_dir, esvee_dir, purple_dir ->
            runnable: bamtools_tumor_dir && purple_dir
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta, redux_tumor_tsv, redux_normal_tsv, bamtools_tumor_dir, bamtools_normal_dir, cobalt_dir, esvee_dir, purple_dir ]
    ch_qsee_inputs = ch_inputs_sorted.runnable
        .map { meta, redux_tumor_tsv, redux_normal_tsv, bamtools_tumor_dir, bamtools_normal_dir, cobalt_dir, esvee_dir, purple_dir ->

            def tumor_id = Utils.getTumorDnaSampleName(meta)
            def normal_id = Utils.getNormalDnaSampleName(meta)

            def meta_qsee = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: tumor_id,
                normal_id: tumor_id,
            ]

            return [ meta, redux_tumor_tsv, redux_normal_tsv, bamtools_tumor_dir, bamtools_normal_dir, cobalt_dir, esvee_dir, purple_dir ]
        }

    // Run process
    QSEE(
        ch_qsee_inputs,
        driver_gene_panel,
        qsee_cohort_percentiles,
    )

    ch_versions = ch_versions.mix(QSEE.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, chord_dir ]
    ch_outputs = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(QSEE.out.qsee_dir, ch_inputs),
            PlaceholderChannels.toolDir(ch_inputs_sorted.skip),
        )

    emit:
    qsee_dir = ch_outputs  // channel: [ meta, qsee_dir ]

    versions  = ch_versions // channel: [ versions.yml ]
}
