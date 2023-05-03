//
// XXX
//
import Constants

include { CUPPA } from '../../modules/local/cuppa/main'

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

        CUPPA(
            ch_cuppa_inputs,
            ref_data_genome_version,
            ref_data_cuppa_resources,
        )

        // Set outputs, restoring original meta
        ch_output = WorkflowOncoanalyser.restoreMeta(CUPPA.out.cuppa_dir, ch_inputs)

    emit:
        cuppa_dir = ch_output

        versions  = CUPPA.out.versions  // channel: [versions.yml]
}
