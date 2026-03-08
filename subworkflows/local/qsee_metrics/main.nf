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

    // Select and route inputs
    // channel: { meta, redux_tumor_tsv, redux_normal_tsv, bamtools_tumor_dir, bamtools_normal_dir, cobalt_dir, esvee_dir, purple_dir }
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
            tumor_bqr_tsv , tumor_dup_freq_tsv , tumor_jitter_tsv , tumor_ms_tsv,
            normal_bqr_tsv, normal_dup_freq_tsv, normal_jitter_tsv, normal_ms_tsv,
            bamtools_tumor_dir, bamtools_normal_dir,
            cobalt_dir,
            esvee_dir,
            purple_dir ->

            def inputs = [:]

            inputs.meta = meta

            inputs.redux_tumor_tsv = [
                Inputs.fallbackToExistingInput(tumor_bqr_tsv, meta, Constants.INPUT.REDUX_BQR_TSV_TUMOR),
                Inputs.fallbackToExistingInput(tumor_dup_freq_tsv, meta, Constants.INPUT.REDUX_DUP_FREQ_TSV_TUMOR),
                Inputs.fallbackToExistingInput(tumor_ms_tsv, meta, Constants.INPUT.REDUX_MS_TSV_TUMOR),
            ].findAll { it != [] }

            inputs.redux_normal_tsv = [
                Inputs.fallbackToExistingInput(normal_bqr_tsv, meta, Constants.INPUT.REDUX_BQR_TSV_NORMAL),
                Inputs.fallbackToExistingInput(normal_dup_freq_tsv, meta, Constants.INPUT.REDUX_DUP_FREQ_TSV_NORMAL),
                Inputs.fallbackToExistingInput(normal_ms_tsv, meta, Constants.INPUT.REDUX_MS_TSV_NORMAL),
            ].findAll { it != [] }

            inputs.bamtools_tumor_dir = Inputs.preferUserProvidedInput(bamtools_tumor_dir, meta, Constants.INPUT.BAMTOOLS_DIR_TUMOR)
            inputs.bamtools_normal_dir = Inputs.preferUserProvidedInput(bamtools_normal_dir, meta, Constants.INPUT.BAMTOOLS_DIR_NORMAL)
            inputs.cobalt_dir = Inputs.preferUserProvidedInput(cobalt_dir, meta, Constants.INPUT.COBALT_DIR)
            inputs.esvee_dir = Inputs.preferUserProvidedInput(esvee_dir, meta, Constants.INPUT.ESVEE_DIR)
            inputs.purple_dir = Inputs.preferUserProvidedInput(purple_dir, meta, Constants.INPUT.PURPLE_DIR)

            return inputs
        }
        .branch { inputs ->
            runnable: inputs.bamtools_tumor_dir && inputs.purple_dir
                return inputs
            skip: true
                return inputs.meta
        }

    // Create process input channel; form metadata
    // channel: [ meta, redux_tumor_tsv, redux_normal_tsv, bamtools_tumor_dir, bamtools_normal_dir, cobalt_dir, esvee_dir, purple_dir ]
    ch_qsee_inputs = ch_inputs_sorted.runnable
        .map { inputs ->

            def meta = inputs.meta

            def meta_qsee = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Inputs.getTumorDnaSampleName(meta),
                normal_id: Inputs.getNormalDnaSampleName(meta),
            ]

            inputs.meta = meta_qsee

            return inputs.values()
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
