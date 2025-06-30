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
    ch_bamtools_somatic         // channel: [mandatory] [ meta, metrics_dir ]
    ch_bamtools_germline        // channel: [mandatory] [ meta, metrics_dir ]
    ch_sage_somatic             // channel: [mandatory] [ meta, sage_dir ]
    ch_sage_germline            // channel: [mandatory] [ meta, sage_dir ]
    ch_sage_somatic_append      // channel: [mandatory] [ meta, sage_append_vcf ]
    ch_sage_germline_append     // channel: [mandatory] [ meta, sage_append_vcf ]
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
    cohort_percentiles          // channel: [mandatory] /path/to/cohort_percentiles
    known_fusion_data           // channel: [mandatory] /path/to/known_fusion_data
    driver_gene_panel           // channel: [mandatory] /path/to/driver_gene_panel
    sigs_etiology               // channel: [mandatory] /path/to/sigs_etiology
    ensembl_data_resources      // channel: [mandatory] /path/to/ensembl_data_resources/
    isofox_alt_sj               // channel: [optional]  /path/to/isofox_alt_sj
    isofox_gene_distribution    // channel: [optional]  /path/to/isofox_gene_distribution

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Set expected input ordering and size
    input_expected_size = 17

    dna_tumor_input_indexes = [
        0,   // bamtools_somatic
        2,   // sage_somatic
        6,   // purple_dir
        7,   // linx_somatic_annotation
        8,   // linx_somatic_plot_dir
    ]

    dna_normal_input_indexes = [
        1,   // bamtools_germline
        3,   // sage_germline
        9,   // linx_germline_annotation
    ]

    rna_tumor_input_indexes = [
        4,   // sage_somatic_append
        16,  // isofox_dir
    ]

    rna_sage_germline_append_index = 7  // sage_germline_append
    cuppa_dir_index = 14                // cuppa_dir

    // Select input sources
    // channel: [ meta, tbt_metrics_dir, nbt_metrics_dir, tsage_dir, nsage_dir, tsage_append, nsage_append, purple_dir, tlinx_anno_dir, tlinx_plot_dir, nlinx_anno_dir, virusinterpreter_dir, chord_dir, sigs_dir, lilac_dir, cuppa_dir, ch_peach, isofox_dir ]
    ch_inputs_selected = WorkflowOncoanalyser.groupByMeta(
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

            assert inputs.size() == input_expected_size

            // NOTE(SW): avoiding further complexity with loops etc

            def inputs_selected = [
                Utils.selectCurrentOrExisting(inputs[0], meta, Constants.INPUT.BAMTOOLS_DIR_TUMOR),
                Utils.selectCurrentOrExisting(inputs[1], meta, Constants.INPUT.BAMTOOLS_DIR_NORMAL),
                Utils.selectCurrentOrExisting(inputs[2], meta, Constants.INPUT.SAGE_DIR_TUMOR),
                Utils.selectCurrentOrExisting(inputs[3], meta, Constants.INPUT.SAGE_DIR_NORMAL),
                Utils.selectCurrentOrExisting(inputs[4], meta, Constants.INPUT.SAGE_APPEND_VCF_TUMOR),
                Utils.selectCurrentOrExisting(inputs[5], meta, Constants.INPUT.SAGE_APPEND_VCF_NORMAL),
                Utils.selectCurrentOrExisting(inputs[6], meta, Constants.INPUT.PURPLE_DIR),
                Utils.selectCurrentOrExisting(inputs[7], meta, Constants.INPUT.LINX_ANNO_DIR_TUMOR),
                Utils.selectCurrentOrExisting(inputs[8], meta, Constants.INPUT.LINX_PLOT_DIR_TUMOR),
                Utils.selectCurrentOrExisting(inputs[9], meta, Constants.INPUT.LINX_ANNO_DIR_NORMAL),
                Utils.selectCurrentOrExisting(inputs[10], meta, Constants.INPUT.VIRUSINTERPRETER_DIR),
                Utils.selectCurrentOrExisting(inputs[11], meta, Constants.INPUT.CHORD_DIR),
                Utils.selectCurrentOrExisting(inputs[12], meta, Constants.INPUT.SIGS_DIR),
                Utils.selectCurrentOrExisting(inputs[13], meta, Constants.INPUT.LILAC_DIR),
                Utils.selectCurrentOrExisting(inputs[14], meta, Constants.INPUT.CUPPA_DIR),
                Utils.selectCurrentOrExisting(inputs[15], meta, Constants.INPUT.PEACH_DIR),
                Utils.selectCurrentOrExisting(inputs[16], meta, Constants.INPUT.ISOFOX_DIR),
            ]

            return [meta, *inputs_selected]
        }

    // Sort inputs
    // channel: runnable: [ meta, tbt_metrics_dir, nbt_metrics_dir, tsage_dir, nsage_dir, tsage_append, nsage_append, purple_dir, tlinx_anno_dir, tlinx_plot_dir, nlinx_anno_dir, virusinterpreter_dir, chord_dir, sigs_dir, lilac_dir, cuppa_dir, peach_dir, isofox_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_inputs_selected
        .branch { d ->

            def meta = d[0]
            def inputs = d[1..-1]

            def has_dna_tumor = dna_tumor_input_indexes
                .collect { i -> inputs[i] }
                .every()

            def has_rna_tumor = rna_tumor_input_indexes
                .collect { i -> inputs[i] }
                .every()

            runnable_dna_and_rna: has_dna_tumor && has_rna_tumor
            runnable_dna: has_dna_tumor
            skip: true
                return meta
        }

    // First set RNA reference files
    // NOTE(SW): since the RNA reference files are provided as channels, I seem to be only able to include via channel ops
    // channel: [ meta, tbt_metrics_dir, nbt_metrics_dir, tsage_dir, nsage_dir, tsage_append, nsage_append, purple_dir, tlinx_anno_dir, tlinx_plot_dir, nlinx_anno_dir, virusinterpreter_dir, chord_dir, sigs_dir, lilac_dir, cuppa_dir, peach_dir, isofox_dir, isofox_alt_sj, isofox_gene_distribution ]
    ch_inputs_runnable = Channel.empty()
        .mix(
            ch_inputs_sorted.runnable_dna.map { d -> [*d, [], []] },
            ch_inputs_sorted.runnable_dna_and_rna
                .combine(isofox_alt_sj)
                .combine(isofox_gene_distribution),
        )

    // Create process input channel
    // channel: sample_data: [ meta, tbt_metrics_dir, nbt_metrics_dir, tsage_dir, nsage_dir, tsmlv_vcf, nsmlv_vcf, purple_dir, tlinx_anno_dir, tlinx_plot_dir, nlinx_anno_dir, virusinterpreter_dir, chord_dir, sigs_dir, lilac_dir, cuppa_dir, peach_dir, isofox_dir ]
    // channel: isofox_alt_sj: [ isofox_alt_sj ]
    // channel: isofox_gene_distribution: [ isofox_gene_distribution ]
    ch_orange_inputs = ch_inputs_runnable
        .multiMap { d ->

            def meta = d[0]
            def inputs = d[1..-3]

            def isofox_alt_sj = d[-2]
            def isofox_gene_distribution = d[-1]

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
                rna_tumor_input_indexes_ready = [*rna_tumor_input_indexes, rna_sage_germline_append_index]
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
            if (inputs_selected[cuppa_dir_index]) {
                def cuppa_vis_data_fp = inputs_selected[cuppa_dir_index].resolve("${meta_orange.tumor_id}.cuppa.vis_data.tsv")
                if (! cuppa_vis_data_fp.exists()) {
                    inputs_selected[cuppa_dir_index] = []
                }
            }

            assert inputs_selected.size() == input_expected_size

            sample_data: [meta_orange, *inputs_selected]
            isofox_alt_sj: isofox_alt_sj
            isofox_gene_distribution: isofox_gene_distribution
        }

    // Run process
    ORANGE(
        ch_orange_inputs.sample_data,
        genome_version,
        disease_ontology,
        cohort_mapping,
        cohort_percentiles,
        known_fusion_data,
        driver_gene_panel,
        sigs_etiology,
        ensembl_data_resources,
        ch_orange_inputs.isofox_alt_sj,
        ch_orange_inputs.isofox_gene_distribution,
        '2.1.0 [oncoanalyser]',
    )

    ch_versions = ch_versions.mix(ORANGE.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
