//
// XXX
//
import Constants

include { CUPPA_CLASSIFIER as CLASSIFIER } from '../../modules/local/cuppa/classifier/main'
include { CUPPA_VISUALISER as VISUALISER } from '../../modules/local/cuppa/visualiser/main'

workflow CUPPA_PREDICTION {
    take:
        // Sample data
        ch_inputs
        ch_inputs_wts_absent
        ch_isofox
        ch_purple
        ch_linx
        ch_virusinterpreter

        // Reference data
        ref_data_genome_version
        ref_data_cuppa_resources

        // Parameters
        run

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [val(meta), isofox_dir]
        ch_cuppa_inputs_isofox = Channel.empty()
        if (run.isofox) {
            // Take Isofox output and supplement with placeholders for missing; i.e. allow optional
            ch_cuppa_inputs_isofox = ch_cuppa_inputs_isofox
                .mix(
                    ch_isofox,
                    ch_inputs_wts_absent.map { meta -> [meta, []] },
                )
        } else {
            ch_cuppa_inputs_isofox = WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.ISOFOX, type: 'optional')
        }

        // channel: [val(meta), isofox_dir, purple_dir, linx_dir, virusinterpreter]
        ch_cuppa_inputs_source = WorkflowOncoanalyser.groupByMeta(
            ch_cuppa_inputs_isofox,
            run.purple ? ch_purple : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PURPLE_DIR, type: 'optional'),
            run.linx ? ch_linx : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.LINX_ANNO_DIR_TUMOR, type: 'optional'),
            run.virusinterpreter ? ch_virusinterpreter : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.VIRUSINTERPRETER_TSV, type: 'optional'),
            flatten_mode: 'nonrecursive',
        )

        // Create inputs and create process-specific meta
        // channel: [val(meta_cuppa), isofox_dir, purple_dir, linx_dir, virusinterpreter]
        ch_cuppa_inputs = ch_cuppa_inputs_source
            .map { data ->
                def meta = data[0]
                def meta_cuppa = [key: meta.id]

                def sample_name_wgs = meta.getAt(['sample_name', Constants.SampleType.TUMOR, Constants.SequenceType.WGS])
                def sample_name_wts = meta.getAt(['sample_name', Constants.SampleType.TUMOR, Constants.SequenceType.WTS])

                if (sample_name_wgs && sample_name_wts) {
                    meta_cuppa.id = sample_name_wgs
                    meta_cuppa.id_wts = sample_name_wts
                } else if (sample_name_wgs) {
                    meta_cuppa.id = sample_name_wgs
                } else if (sample_name_wts) {
                    meta_cuppa.id = sample_name_wts
                } else {
                    log.error "ERROR: no sample name for: ${meta}"
                    System.exit(1)
                }

                return [meta_cuppa, *data[1..-1]]
            }

        CLASSIFIER(
            ch_cuppa_inputs,
            ref_data_genome_version,
            ref_data_cuppa_resources,
        )

        VISUALISER(
            CLASSIFIER.out.csv,
        )

        // Set outputs, restoring original meta
        ch_csv = WorkflowOncoanalyser.restoreMeta(CLASSIFIER.out.csv, ch_inputs)
        ch_summary_plot = WorkflowOncoanalyser.restoreMeta(VISUALISER.out.summary_plot, ch_inputs)
        ch_feature_plot = WorkflowOncoanalyser.restoreMeta(VISUALISER.out.feature_plot, ch_inputs)
        ch_versions = ch_versions.mix(
            VISUALISER.out.versions,
            CLASSIFIER.out.versions,
        )

    emit:
        csv = ch_csv
        summary_plot = ch_summary_plot
        feature_plot = ch_feature_plot

        versions  = ch_versions // channel: [versions.yml]
}
