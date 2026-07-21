//
// ORANGE collates outputs of hmftools into a static PDF report
//

import Constants
import Utils

include { ORANGE } from '../../../modules/local/orange/main'

workflow ORANGE_REPORTING {
    take:
    // Sample data
    ch_sage_dir_somatic             // channel: [mandatory] [ meta, sage_dir ]
    ch_sage_dir_germline            // channel: [mandatory] [ meta, sage_dir ]
    ch_sage_append_dir_somatic      // channel: [mandatory] [ meta, sage_append_dir ]
    ch_sage_append_dir_germline     // channel: [mandatory] [ meta, sage_append_dir ]
    ch_sage_plot_dir_somatic        // channel: [mandatory] [ meta, sage_visualiser_dir ]
    ch_purple_dir                   // channel: [mandatory] [ meta, purple_dir ]
    ch_qsee_dir                     // channel: [mandatory] [ meta, qsee_dir ]
    ch_linx_annotation_dir_somatic  // channel: [mandatory] [ meta, linx_annotation_dir ]
    ch_linx_plot_dir_somatic        // channel: [mandatory] [ meta, linx_visualiser_dir ]
    ch_linx_annotation_dir_germline // channel: [mandatory] [ meta, linx_annotation_dir ]
    ch_virusinterpreter_dir         // channel: [mandatory] [ meta, virusinterpreter_dir ]
    ch_chord_dir                    // channel: [mandatory] [ meta, chord_dir ]
    ch_sigs_dir                     // channel: [mandatory] [ meta, sigs_dir ]
    ch_lilac_dir                    // channel: [optional]  [ meta, lilac_dir ]
    ch_cuppa_dir                    // channel: [mandatory] [ meta, cuppa_dir ]
    ch_peach_dir                    // channel: [mandatory] [ meta, peach_dir ]
    ch_isofox_dir                   // channel: [mandatory] [ meta, isofox_dir ]

    // Reference data
    genome_version                  // channel: [mandatory] genome version
    disease_ontology                // channel: [mandatory] /path/to/disease_ontology

    // Params
    sequencing_platform             // string:  [mandatory] sequencing platform
    targeted_mode                   // boolean: [mandatory] Set targeted mode

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Mapping for semantic input retrieval
    input_indexes = [
        'sage_dir_somatic':             0,
        'sage_dir_germline':            1,
        'sage_append_dir_somatic':      2,
        'sage_append_dir_germline':     3,
        'sage_plot_dir_somatic':        4,
        'purple_dir':                   5,
        'qsee_dir':                     6,
        'linx_annotation_dir_somatic':  7,
        'linx_plot_dir_somatic':        8,
        'linx_annotation_dir_germline': 9,
        'virusinterpreter_dir':         10,
        'chord_dir':                    11,
        'sigs_dir':                     12,
        'lilac_dir':                    13,
        'cuppa_dir':                    14,
        'peach_dir':                    15,
        'isofox_dir':                   16,
    ]

    // Select input sources then sort
    // channel: runnable: [meta, sage_dir_somatic, sage_dir_germline, sage_append_dir_somatic, sage_append_dir_germline, sage_plot_dir_somatic, purple_dir, qsee_dir, linx_annotation_dir_somatic, linx_plot_dir_somatic, linx_annotation_dir_germline, virusinterpreter_dir, chord_dir, sigs_dir, lilac_dr, cuppa_dir, peach_dir, isofox_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_sage_dir_somatic,
        ch_sage_dir_germline,
        ch_sage_append_dir_somatic,
        ch_sage_append_dir_germline,
        ch_sage_plot_dir_somatic,
        ch_purple_dir,
        ch_qsee_dir,
        ch_linx_annotation_dir_somatic,
        ch_linx_plot_dir_somatic,
        ch_linx_annotation_dir_germline,
        ch_virusinterpreter_dir,
        ch_chord_dir,
        ch_sigs_dir,
        ch_lilac_dir,
        ch_cuppa_dir,
        ch_peach_dir,
        ch_isofox_dir,
    )
        .map { d ->

            def meta = d[0]
            def inputs = d[1..-1]

            assert inputs.size() == input_indexes.size()

            // NOTE(SW): avoiding further complexity with loops etc
            return [
                meta,
                Utils.selectCurrentOrExisting(inputs[input_indexes['sage_dir_somatic']],             meta, Constants.INPUT.SAGE_DIR_TUMOR),
                Utils.selectCurrentOrExisting(inputs[input_indexes['sage_dir_germline']],            meta, Constants.INPUT.SAGE_DIR_NORMAL),
                Utils.selectCurrentOrExisting(inputs[input_indexes['sage_append_dir_somatic']],      meta, Constants.INPUT.SAGE_APPEND_DIR_TUMOR),
                Utils.selectCurrentOrExisting(inputs[input_indexes['sage_append_dir_germline']],     meta, Constants.INPUT.SAGE_APPEND_DIR_NORMAL),
                Utils.selectCurrentOrExisting(inputs[input_indexes['sage_plot_dir_somatic']],        meta, Constants.INPUT.SAGE_PLOT_DIR_TUMOR),
                Utils.selectCurrentOrExisting(inputs[input_indexes['purple_dir']],                   meta, Constants.INPUT.PURPLE_DIR),
                Utils.selectCurrentOrExisting(inputs[input_indexes['qsee_dir']],                     meta, Constants.INPUT.QSEE_DIR),
                Utils.selectCurrentOrExisting(inputs[input_indexes['linx_annotation_dir_somatic']],  meta, Constants.INPUT.LINX_ANNO_DIR_TUMOR),
                Utils.selectCurrentOrExisting(inputs[input_indexes['linx_plot_dir_somatic']],        meta, Constants.INPUT.LINX_PLOT_DIR_TUMOR),
                Utils.selectCurrentOrExisting(inputs[input_indexes['linx_annotation_dir_germline']], meta, Constants.INPUT.LINX_ANNO_DIR_NORMAL),
                Utils.selectCurrentOrExisting(inputs[input_indexes['virusinterpreter_dir']],         meta, Constants.INPUT.VIRUSINTERPRETER_DIR),
                Utils.selectCurrentOrExisting(inputs[input_indexes['chord_dir']],                    meta, Constants.INPUT.CHORD_DIR),
                Utils.selectCurrentOrExisting(inputs[input_indexes['sigs_dir']],                     meta, Constants.INPUT.SIGS_DIR),
                Utils.selectCurrentOrExisting(inputs[input_indexes['lilac_dir']],                    meta, Constants.INPUT.LILAC_DIR),
                Utils.selectCurrentOrExisting(inputs[input_indexes['cuppa_dir']],                    meta, Constants.INPUT.CUPPA_DIR),
                Utils.selectCurrentOrExisting(inputs[input_indexes['peach_dir']],                    meta, Constants.INPUT.PEACH_DIR),
                Utils.selectCurrentOrExisting(inputs[input_indexes['isofox_dir']],                   meta, Constants.INPUT.ISOFOX_DIR),
            ]

        }
        .branch { d ->

            def meta = d[0]
            def inputs = d[1..-1]

            def dna_tumor_input_keys = ['sage_dir_somatic', 'purple_dir', 'qsee_dir', 'linx_annotation_dir_somatic', 'linx_plot_dir_somatic']
            def has_dna_tumor = dna_tumor_input_keys.every { k -> i = input_indexes[k]; return inputs[i] }

            runnable: has_dna_tumor
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [meta, sage_dir_somatic, sage_dir_germline, sage_append_dir_somatic, sage_append_dir_germline, sage_plot_dir_somatic, purple_dir, qsee_dir, linx_annotation_dir_somatic, linx_plot_dir_somatic, linx_annotation_dir_germline, virusinterpreter_dir, chord_dir, sigs_dir, lilac_dr, cuppa_dir, peach_dir, isofox_dir ]
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
            def dna_normal_input_keys = ['sage_dir_germline', 'linx_annotation_dir_germline']
            def has_dna_normal = dna_normal_input_keys.every { k -> i = input_indexes[k]; return inputs[i] }

            if (has_dna_normal) {
                meta_orange.normal_dna_id = Utils.getNormalDnaSampleName(meta)
            } else {
                dna_normal_input_keys.each { k -> i = input_indexes[k]; inputs_selected[i] = [] }
            }

            // Require all tumor RNA inputs to be present else clear them
            // SAGE append germline is only required when normal DNA is present
            def rna_tumor_input_keys
            if (has_dna_normal) {
                rna_tumor_input_keys = ['isofox_dir', 'sage_append_dir_somatic', 'sage_append_dir_germline']
            } else {
                rna_tumor_input_keys = ['isofox_dir', 'sage_append_dir_somatic']
            }

            def has_rna_tumor = rna_tumor_input_keys.every { k -> i = input_indexes[k]; return inputs[i] }

            if (has_rna_tumor) {
                meta_orange.tumor_rna_id = Utils.getTumorRnaSampleName(meta)
            } else {
                rna_tumor_input_keys.each { k -> i = input_indexes[k]; inputs_selected[i] = [] }
            }

            // ORANGE only accepts CUPPA with DNA; when providing DNA/RNA inputs but skipping Virus Interpreter CUPPA
            // will generate RNA only outputs and no visualisation, which triggers missing file error in ORANGE
            def cuppa_dir = inputs[input_indexes['cuppa_dir']]
            if (cuppa_dir) {
                def cuppa_vis_data = cuppa_dir.resolve("${meta_orange.tumor_id}.cuppa.vis_data.tsv")
                if (! cuppa_vis_data.exists()) {
                    inputs_selected[input_indexes['cuppa_dir']] = []
                }
            }

            // Set SAGE append VCF input
            def sage_append_vcf_somatic
            def sage_append_vcf_germline
            if (has_rna_tumor) {
                // Somatic
                def sage_append_dir_somatic = inputs_selected[input_indexes['sage_append_dir_somatic']]
                if (sage_append_dir_somatic) {
                    sage_append_vcf_somatic = sage_append_dir_somatic.resolve("${meta_orange.tumor_id}.sage.append.vcf.gz")
                }

                // Germline
                def sage_append_dir_germline = inputs_selected[input_indexes['sage_append_dir_germline']]
                if (sage_append_dir_germline) {
                    sage_append_vcf_germline = sage_append_dir_germline.resolve("${meta_orange.normal_dna_id}.sage.append.vcf.gz")
                }
            }

            // Set LINX reportable plot directory
            def linx_plot_dir_somatic = inputs_selected[input_indexes['linx_plot_dir_somatic']]
            def linx_plot_dir_somatic_reportable = linx_plot_dir_somatic.resolve('reportable/')

            // The LINX directory may not exist on object store providers where no plots where created
            linx_plot_dir_somatic_reportable = linx_plot_dir_somatic_reportable.exists() ? linx_plot_dir_somatic_reportable : []

            return [
                meta_orange,
                inputs_selected[input_indexes['sage_dir_somatic']],
                inputs_selected[input_indexes['sage_dir_germline']],
                sage_append_vcf_somatic ?: [],
                sage_append_vcf_germline ?: [],
                inputs_selected[input_indexes['sage_plot_dir_somatic']],
                inputs_selected[input_indexes['purple_dir']],
                inputs_selected[input_indexes['qsee_dir']],
                inputs_selected[input_indexes['linx_annotation_dir_somatic']],
                linx_plot_dir_somatic_reportable,
                inputs_selected[input_indexes['linx_annotation_dir_germline']],
                inputs_selected[input_indexes['virusinterpreter_dir']],
                inputs_selected[input_indexes['chord_dir']],
                inputs_selected[input_indexes['sigs_dir']],
                inputs_selected[input_indexes['lilac_dir']],
                inputs_selected[input_indexes['cuppa_dir']],
                inputs_selected[input_indexes['peach_dir']],
                inputs_selected[input_indexes['isofox_dir']],
            ]
        }

    // Run process
    ORANGE(
        ch_orange_inputs,
        genome_version,
        disease_ontology,
        '2.3.0 [oncoanalyser]',
<<<<<<< HEAD
        sequencing_type,
=======
        sequencing_platform,
>>>>>>> ed465085 (Hartwig integration with backported fixes)
        targeted_mode,
    )

    ch_versions = ch_versions.mix(ORANGE.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
