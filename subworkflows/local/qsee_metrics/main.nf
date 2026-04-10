//
// Qsee calculates and visualises QC metrics
//

include { QSEE } from '../../../modules/local/qsee/main'

workflow QSEE_METRICS {
    take:
    // Sample data
    ch_inputs                // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor       // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal      // channel: [mandatory] [ meta, redux_dir ]
    ch_bamtools_dir_tumor    // channel: [mandatory] [ meta, metrics_dir ]
    ch_bamtools_dir_normal   // channel: [optional]  [ meta, metrics_dir ]
    ch_cobalt_dir            // channel: [optional]  [ meta, cobalt_dir ]
    ch_esvee_dir             // channel: [optional]  [ meta, esvee_dir ]
    ch_purple_dir            // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    driver_gene_panel        // channel: [mandatory] /path/to/driver_gene_panel
    qsee_cohort_percentiles  // channel: [mandatory] /path/to/cohort_percentiles

    // Params
    sequencing_type          // string:  [mandatory] sequencing type
    targeted_mode            // boolean: [mandatory] Set targeted mode

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select and route inputs
    // channel: { meta, redux_dir_tumor, redux_dir_normal, bamtools_tumor_dir, bamtools_normal_dir, cobalt_dir, esvee_dir, purple_dir }
    ch_inputs_sorted = channels.WorkflowChannels.groupByMeta(
        flatten_mode: 'singletons_only',
        ch_redux_dir_tumor,
        ch_redux_dir_normal,
        ch_bamtools_dir_tumor,
        ch_bamtools_dir_normal,
        ch_cobalt_dir,
        ch_esvee_dir,
        ch_purple_dir,
    )
        .map { meta,
            redux_tsvs_tumor, redux_tsvs_normal,
            bamtools_tumor_dir, bamtools_normal_dir,
            cobalt_dir,
            esvee_dir,
            purple_dir ->

            def inputs = [:]

            inputs.meta = meta

            inputs.redux_tsvs_tumor = sample.Inputs.resolveReduxTsvFiles(redux_tsvs_tumor, meta, samplesheet.SampleType.TUMOR)
            inputs.redux_tsvs_normal = sample.Inputs.resolveReduxTsvFiles(redux_tsvs_normal, meta, samplesheet.SampleType.NORMAL)

            inputs.bamtools_tumor_dir = sample.Inputs.preferUserProvidedInput(bamtools_tumor_dir, meta, sample.FileKey.BAMTOOLS_DIR_TUMOR)
            inputs.bamtools_normal_dir = sample.Inputs.preferUserProvidedInput(bamtools_normal_dir, meta, sample.FileKey.BAMTOOLS_DIR_NORMAL)

            inputs.cobalt_dir = sample.Inputs.preferUserProvidedInput(cobalt_dir, meta, sample.FileKey.COBALT_DIR)
            inputs.esvee_dir = sample.Inputs.preferUserProvidedInput(esvee_dir, meta, sample.FileKey.ESVEE_DIR)
            inputs.purple_dir = sample.Inputs.preferUserProvidedInput(purple_dir, meta, sample.FileKey.PURPLE_DIR)

            return inputs
        }
        .branch { inputs ->
            runnable: inputs.bamtools_tumor_dir && inputs.purple_dir
                return inputs
            skip: true
                return inputs.meta
        }

    // Create process input channel; form metadata
    // channel: [ meta, redux_tsvs_tumor, redux_tsvs_normal, bamtools_tumor_dir, bamtools_normal_dir, cobalt_dir, esvee_dir, purple_dir ]
    ch_qsee_inputs = ch_inputs_sorted.runnable
        .map { inputs ->

            def meta = inputs.meta

            def has_normal_dna = inputs.redux_tsvs_normal || inputs.bamtools_normal_dir

            def meta_qsee = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: sample.Inputs.getTumorDnaSampleName(meta),
                normal_id: has_normal_dna ? sample.Inputs.getNormalDnaSampleName(meta) : null,
            ]

            inputs.meta = meta_qsee

            return inputs.values()
        }

    // Run process
    QSEE(
        ch_qsee_inputs,
        driver_gene_panel,
        qsee_cohort_percentiles,
        sequencing_type,
        targeted_mode,
    )

    ch_versions = ch_versions.mix(QSEE.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, chord_dir ]
    ch_outputs = Channel.empty()
        .mix(
            channels.WorkflowChannels.restoreMeta(QSEE.out.qsee_dir, ch_inputs),
            channels.PlaceholderChannels.toolDir(ch_inputs_sorted.skip),
        )

    emit:
    qsee_dir = ch_outputs  // channel: [ meta, qsee_dir ]

    versions  = ch_versions // channel: [ versions.yml ]
}
