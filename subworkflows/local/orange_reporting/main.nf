//
// ORANGE collates outputs of hmftools into a static PDF report
//

import Constants
import Utils

include { ORANGE } from '../../../modules/local/orange/main'

workflow ORANGE_REPORTING {
    take:
    // Sample data
    ch_inputs                   // channel: [mandatory] [ meta ]
    ch_redux_somatic_plot       // channel: [mandatory] [ meta, redux_bqr_plot ]
    ch_redux_germline_plot      // channel: [mandatory] [ meta, redux_bqr_plot ]
    ch_bamtools_somatic         // channel: [mandatory] [ meta, metrics_dir ]
    ch_bamtools_germline        // channel: [mandatory] [ meta, metrics_dir ]
    ch_sage_somatic             // channel: [mandatory] [ meta, sage_dir ]
    ch_sage_germline            // channel: [mandatory] [ meta, sage_dir ]
    ch_sage_somatic_append      // channel: [mandatory] [ meta, sage_append_dir ]
    ch_sage_germline_append     // channel: [mandatory] [ meta, sage_append_dir ]
    ch_purple                   // channel: [mandatory] [ meta, purple_dir ]
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
    cohort_mapping              // channel: [mandatory] /path/to/cohort_mapping
    driver_gene_panel           // channel: [mandatory] /path/to/driver_gene_panel
    sigs_etiology               // channel: [mandatory] /path/to/sigs_etiology
    ensembl_data_resources      // channel: [mandatory] /path/to/ensembl_data_resources/

    // Params
    targeted_mode               // boolean: [mandatory] Set targeted mode

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Refer to inputs by index to avoid needing to declare each variable when calling map, branch, etc
    input_indexes = [
        'redux_somatic_plot'      : 0,
        'redux_germline_plot'     : 1,
        'bamtools_somatic'        : 2,
        'bamtools_germline'       : 3,
        'sage_somatic'            : 4,
        'sage_germline'           : 5,
        'sage_somatic_append'     : 6,
        'sage_germline_append'    : 7,
        'purple_dir'              : 8,
        'linx_somatic_annotation' : 9,
        'linx_somatic_plot'       : 10,
        'linx_germline_annotation': 11,
        'virusinterpreter'        : 12,
        'chord'                   : 13,
        'sigs'                    : 14,
        'lilac'                   : 15,
        'cuppa'                   : 16,
        'peach'                   : 17,
        'isofox'                  : 18,
    ]

    dna_tumor_input_indexes = [
        input_indexes.redux_somatic_plot,
        input_indexes.bamtools_somatic,
        input_indexes.sage_somatic,
        input_indexes.purple_dir,
        input_indexes.linx_somatic_annotation,
        input_indexes.linx_somatic_plot,
    ]

    dna_normal_input_indexes = [
        input_indexes.redux_germline_plot,
        input_indexes.bamtools_germline,
        input_indexes.sage_germline,
        input_indexes.linx_germline_annotation,
    ]

    rna_tumor_input_indexes = [
        input_indexes.sage_somatic_append,
        input_indexes.isofox,
    ]

    // Select input sources
    // channel: [ meta, redux_somatic_plot, redux_germline_plot, ... ]
    ch_inputs_selected = WorkflowOncoanalyser.groupByMeta(
        ch_redux_somatic_plot,
        ch_redux_germline_plot,
        ch_bamtools_somatic,
        ch_bamtools_germline,
        ch_sage_somatic,
        ch_sage_germline,
        ch_sage_somatic_append,
        ch_sage_germline_append,
        ch_purple,
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
        .map { d ->

            def meta = d[0]
            def inputs = d[1..-1]

            assert inputs.size() == input_indexes.size()

            def inputs_selected = [
                Utils.selectCurrentOrExisting(inputs[input_indexes.redux_somatic_plot]      , meta, Constants.INPUT.REDUX_BQR_PLOT_TUMOR),
                Utils.selectCurrentOrExisting(inputs[input_indexes.redux_germline_plot]     , meta, Constants.INPUT.REDUX_BQR_PLOT_NORMAL),
                Utils.selectCurrentOrExisting(inputs[input_indexes.bamtools_somatic]        , meta, Constants.INPUT.BAMTOOLS_DIR_TUMOR),
                Utils.selectCurrentOrExisting(inputs[input_indexes.bamtools_germline]       , meta, Constants.INPUT.BAMTOOLS_DIR_NORMAL),
                Utils.selectCurrentOrExisting(inputs[input_indexes.sage_somatic]            , meta, Constants.INPUT.SAGE_DIR_TUMOR),
                Utils.selectCurrentOrExisting(inputs[input_indexes.sage_germline]           , meta, Constants.INPUT.SAGE_DIR_NORMAL),
                Utils.selectCurrentOrExisting(inputs[input_indexes.sage_somatic_append]     , meta, Constants.INPUT.SAGE_APPEND_DIR_TUMOR),
                Utils.selectCurrentOrExisting(inputs[input_indexes.sage_germline_append]    , meta, Constants.INPUT.SAGE_APPEND_DIR_NORMAL),
                Utils.selectCurrentOrExisting(inputs[input_indexes.purple_dir]              , meta, Constants.INPUT.PURPLE_DIR),
                Utils.selectCurrentOrExisting(inputs[input_indexes.linx_somatic_annotation] , meta, Constants.INPUT.LINX_ANNO_DIR_TUMOR),
                Utils.selectCurrentOrExisting(inputs[input_indexes.linx_somatic_plot]       , meta, Constants.INPUT.LINX_PLOT_DIR_TUMOR),
                Utils.selectCurrentOrExisting(inputs[input_indexes.linx_germline_annotation], meta, Constants.INPUT.LINX_ANNO_DIR_NORMAL),
                Utils.selectCurrentOrExisting(inputs[input_indexes.virusinterpreter]        , meta, Constants.INPUT.VIRUSINTERPRETER_DIR),
                Utils.selectCurrentOrExisting(inputs[input_indexes.chord]                   , meta, Constants.INPUT.CHORD_DIR),
                Utils.selectCurrentOrExisting(inputs[input_indexes.sigs]                    , meta, Constants.INPUT.SIGS_DIR),
                Utils.selectCurrentOrExisting(inputs[input_indexes.lilac]                   , meta, Constants.INPUT.LILAC_DIR),
                Utils.selectCurrentOrExisting(inputs[input_indexes.cuppa]                   , meta, Constants.INPUT.CUPPA_DIR),
                Utils.selectCurrentOrExisting(inputs[input_indexes.peach]                   , meta, Constants.INPUT.PEACH_DIR),
                Utils.selectCurrentOrExisting(inputs[input_indexes.isofox]                  , meta, Constants.INPUT.ISOFOX_DIR),
            ]

            return [meta, *inputs_selected]
        }

    // Sort inputs
    // channel: runnable: [ meta, redux_somatic_plot, redux_germline_plot, ... ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_inputs_selected
        .branch { d ->

            def meta = d[0]
            def inputs = d[1..-1]

            def has_dna_tumor = dna_tumor_input_indexes
                .collect { i -> inputs[i] }
                .every()

            runnable: has_dna_tumor
            skip: true
                return meta
        }

    // Create process input channel
    // channel: sample_data: [ meta, redux_somatic_plot, redux_germline_plot, ... ]
    ch_orange_inputs = ch_inputs_sorted.runnable
        .map { d ->

            def meta = d[0]
            def inputs = d[1..-1]

            def meta_orange = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Utils.getTumorDnaSampleName(meta),
                cancer_type: meta[Constants.InfoField.CANCER_TYPE],
            ]

            def inputs_selected = inputs.clone()

            // Require all normal DNA inputs to be present else clear them
            def has_dna_normal = dna_normal_input_indexes
                .collect { i -> inputs[i] }
                .every()

            if (has_dna_normal) {
                meta_orange.normal_dna_id = Utils.getNormalDnaSampleName(meta)
            } else {
                dna_normal_input_indexes.each { i -> inputs_selected[i] = [] }
            }

            // Require all tumor RNA inputs to be present else clear them
            // SAGE append germline is only required when normal DNA is present
            def rna_tumor_input_indexes_ready
            if (has_dna_normal) {
                rna_tumor_input_indexes_ready = [*rna_tumor_input_indexes, input_indexes.sage_germline_append]
            } else {
                rna_tumor_input_indexes_ready = rna_tumor_input_indexes.clone()
            }

            def has_rna_tumor = rna_tumor_input_indexes_ready
                .collect { i -> inputs[i] }
                .every()

            if (has_rna_tumor) {
                meta_orange.tumor_rna_id = Utils.getTumorRnaSampleName(meta)
            } else {
                rna_tumor_input_indexes.each { i -> inputs_selected[i] = [] }
            }

            // ORANGE only accepts CUPPA with DNA; when providing DNA/RNA inputs but skipping Virus Interpreter CUPPA
            // will generate RNA only outputs and no visualisation, which triggers missing file error in ORANGE
            if (inputs_selected[input_indexes.cuppa]) {
                def cuppa_vis_data_fp = inputs_selected[input_indexes.cuppa].resolve("${meta_orange.tumor_id}.cuppa.vis_data.tsv")
                if (!cuppa_vis_data_fp.exists()) {
                    inputs_selected[input_indexes.cuppa] = []
                }
            }

            // Set SAGE append VCF input
            if (has_rna_tumor) {
                // Somatic
                def sage_somatic_append = inputs_selected[input_indexes.sage_somatic_append]
                if (sage_somatic_append) {
                    inputs_selected[input_indexes.sage_somatic_append] = file(sage_somatic_append).resolve("${meta_orange.tumor_id}.sage.append.vcf.gz")
                }

                // Germline
                def sage_germline_append = inputs_selected[input_indexes.sage_germline_append]
                if (sage_germline_append) {
                    inputs_selected[input_indexes.sage_germline_append] = file(sage_germline_append).resolve("${meta_orange.normal_dna_id}.sage.append.vcf.gz")
                }
            }

            assert inputs_selected.size() == input_indexes.size()

            return [meta_orange, *inputs_selected]
        }

    // Run process
    ORANGE(
        ch_orange_inputs,
        genome_version,
        disease_ontology,
        cohort_mapping,
        driver_gene_panel,
        sigs_etiology,
        ensembl_data_resources,
        '3.0.0 [oncoanalyser]',
        targeted_mode,
    )

    ch_versions = ch_versions.mix(ORANGE.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
