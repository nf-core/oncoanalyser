//
// ORANGE collates outputs of hmftools into a static PDF report
//

include { ORANGE } from '../../../modules/local/orange/main'

workflow ORANGE_REPORTING {
    take:
    // Sample data
    ch_sage_somatic             // channel: [mandatory] [ meta, sage_dir ]
    ch_sage_germline            // channel: [mandatory] [ meta, sage_dir ]
    ch_sage_somatic_append      // channel: [mandatory] [ meta, sage_append_dir ]
    ch_sage_germline_append     // channel: [mandatory] [ meta, sage_append_dir ]
    ch_purple                   // channel: [mandatory] [ meta, purple_dir ]
    ch_qsee                     // channel: [mandatory] [ meta, qsee_dir ]
    ch_linx_somatic_annotation  // channel: [mandatory] [ meta, linx_annotation_dir ]
    ch_linx_somatic_plot        // channel: [mandatory] [ meta, linx_visualiser_dir ]
    ch_linx_germline_annotation // channel: [mandatory] [ meta, linx_annotation_dir ]
    ch_virusinterpreter         // channel: [mandatory] [ meta, virusinterpreter_dir ]
    ch_chord                    // channel: [mandatory] [ meta, chord_dir ]
    ch_sigs                     // channel: [mandatory] [ meta, sigs_dir ]
    ch_lilac                    // channel: [optional]  [ meta, lilac_dir ]
    ch_cuppa                    // channel: [mandatory] [ meta, cuppa_dir ]
    ch_peach                    // channel: [mandatory] [ meta, peach_dir ]
    ch_isofox                   // channel: [mandatory] [ meta, isofox_dir ]

    // Reference data
    genome_version              // channel: [mandatory] genome version
    disease_ontology            // channel: [mandatory] /path/to/disease_ontology

    // Params
    targeted_mode               // boolean: [mandatory] Set targeted mode

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources
    // channel: { meta, ... }
    ch_inputs_selected = WorkflowChannels.groupByMeta(
        ch_sage_somatic,
        ch_sage_germline,
        ch_sage_somatic_append,
        ch_sage_germline_append,
        ch_purple,
        ch_qsee,
        ch_linx_somatic_annotation,
        ch_linx_somatic_plot,
        ch_linx_germline_annotation,
        ch_virusinterpreter,
        ch_chord,
        ch_sigs,
        ch_lilac,
        ch_cuppa,
        ch_peach,
        ch_isofox,
    )
        .map { inputs_list ->

            def inputs_map = [:]

            def meta = inputs_list[0]
            inputs_map.meta = meta

            inputs_map.sage_somatic             = Inputs.preferUserProvidedInput(inputs_list[1], meta, SampleSheetFields.INPUT.SAGE_DIR_TUMOR)
            inputs_map.sage_germline            = Inputs.preferUserProvidedInput(inputs_list[2], meta, SampleSheetFields.INPUT.SAGE_DIR_NORMAL)
            inputs_map.sage_somatic_append      = Inputs.preferUserProvidedInput(inputs_list[3], meta, SampleSheetFields.INPUT.SAGE_APPEND_DIR_TUMOR)
            inputs_map.sage_germline_append     = Inputs.preferUserProvidedInput(inputs_list[4], meta, SampleSheetFields.INPUT.SAGE_APPEND_DIR_NORMAL)
            inputs_map.purple_dir               = Inputs.preferUserProvidedInput(inputs_list[5], meta, SampleSheetFields.INPUT.PURPLE_DIR)
            inputs_map.qsee_dir                 = Inputs.preferUserProvidedInput(inputs_list[6], meta, SampleSheetFields.INPUT.QSEE_DIR)
            inputs_map.linx_somatic_annotation  = Inputs.preferUserProvidedInput(inputs_list[7], meta, SampleSheetFields.INPUT.LINX_ANNO_DIR_TUMOR)
            inputs_map.linx_somatic_plot        = Inputs.preferUserProvidedInput(inputs_list[8], meta, SampleSheetFields.INPUT.LINX_PLOT_DIR_TUMOR)
            inputs_map.linx_germline_annotation = Inputs.preferUserProvidedInput(inputs_list[9], meta, SampleSheetFields.INPUT.LINX_ANNO_DIR_NORMAL)
            inputs_map.virusinterpreter         = Inputs.preferUserProvidedInput(inputs_list[10], meta, SampleSheetFields.INPUT.VIRUSINTERPRETER_DIR)
            inputs_map.chord                    = Inputs.preferUserProvidedInput(inputs_list[11], meta, SampleSheetFields.INPUT.CHORD_DIR)
            inputs_map.sigs                     = Inputs.preferUserProvidedInput(inputs_list[12], meta, SampleSheetFields.INPUT.SIGS_DIR)
            inputs_map.lilac                    = Inputs.preferUserProvidedInput(inputs_list[13], meta, SampleSheetFields.INPUT.LILAC_DIR)
            inputs_map.cuppa                    = Inputs.preferUserProvidedInput(inputs_list[14], meta, SampleSheetFields.INPUT.CUPPA_DIR)
            inputs_map.peach                    = Inputs.preferUserProvidedInput(inputs_list[15], meta, SampleSheetFields.INPUT.PEACH_DIR)
            inputs_map.isofox                   = Inputs.preferUserProvidedInput(inputs_list[16], meta, SampleSheetFields.INPUT.ISOFOX_DIR)

            return inputs_map
        }

    def has_required_inputs = { inputs, keys ->
        return keys.every { key -> inputs[key] }
    }

    def clear_inputs = { inputs, keys ->
        keys.each { key -> inputs[key] = [] }
    }

    // Sort inputs
    // channel: runnable: { meta, ... }
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_inputs_selected
        .branch { inputs ->

            def dna_tumor_input_keys = ['sage_somatic', 'purple_dir', 'qsee_dir', 'linx_somatic_annotation', 'linx_somatic_plot']
            def has_dna_tumor = has_required_inputs(inputs, dna_tumor_input_keys)

            runnable: has_dna_tumor
                return inputs
            skip: true
                return inputs.meta
        }

    // Create process input channel
    // channel: sample_data: [ meta, ... ]
    ch_orange_inputs = ch_inputs_sorted.runnable
        .map { inputs ->

            def meta = inputs.meta

            def meta_orange = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Inputs.getTumorDnaSampleName(meta),
                cancer_type: meta[SampleSheetFields.InfoField.CANCER_TYPE],
            ]

            // Require all normal DNA inputs to be present else clear them
            def dna_normal_input_keys = ['sage_germline', 'linx_germline_annotation' ]
            def has_dna_normal = has_required_inputs(inputs, dna_normal_input_keys)

            if (has_dna_normal) {
                meta_orange.normal_dna_id = Inputs.getNormalDnaSampleName(meta)
            } else {
                clear_inputs(inputs, dna_normal_input_keys)
            }

            // ORANGE only accepts CUPPA with DNA; when providing DNA/RNA inputs but skipping Virus Interpreter CUPPA
            // will generate RNA only outputs and no visualisation, which triggers missing file error in ORANGE
            if (inputs.cuppa) {
                def cuppa_vis_data = inputs.cuppa.resolve("${meta_orange.tumor_id}.cuppa.vis_data.tsv")
                if (!cuppa_vis_data.exists()) {
                    inputs.cuppa = []
                }
            }

            // Require all tumor RNA inputs to be present else clear them
            // SAGE append germline is only required when normal DNA is present
            def rna_tumor_input_keys = has_dna_normal
                ? ['isofox', 'sage_somatic_append', 'sage_germline_append']
                : ['isofox', 'sage_somatic_append']

            def has_rna_tumor = has_required_inputs(inputs, rna_tumor_input_keys)

            if (has_rna_tumor) {
                meta_orange.tumor_rna_id = Inputs.getTumorRnaSampleName(meta)
            } else {
                clear_inputs(inputs, rna_tumor_input_keys)
            }

            // Set SAGE append VCF input
            if (has_rna_tumor) {
                def sage_somatic_append_dir = inputs.sage_somatic_append
                if (sage_somatic_append_dir) {
                    inputs.sage_somatic_append = file(sage_somatic_append_dir).resolve("${meta_orange.tumor_id}.sage.append.vcf.gz")
                }

                def sage_germline_append_dir = inputs.sage_germline_append
                if (sage_germline_append_dir) {
                    inputs.sage_germline_append = file(sage_germline_append_dir).resolve("${meta_orange.normal_dna_id}.sage.append.vcf.gz")
                }
            }

            inputs.meta = meta_orange

            return inputs.values()
        }

    // Run process
    ORANGE(
        ch_orange_inputs,
        genome_version,
        disease_ontology,
        '3.0.0 [oncoanalyser]',
        targeted_mode,
    )

    ch_versions = ch_versions.mix(ORANGE.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
