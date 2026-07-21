//
// SAGE append adds additional sample data to an existing SAGE VCF
//

import Constants
import Utils

include { SAGE_APPEND as SAGE_APPEND_SOMATIC  } from '../../../modules/local/sage/append/main'
include { SAGE_APPEND as SAGE_APPEND_GERMLINE } from '../../../modules/local/sage/append/main'

workflow SAGE_APPEND {
    take:
    // Sample data
    ch_inputs           // channel: [mandatory] [ meta ]
    ch_purple_dir       // channel: [mandatory] [ meta, purple_dir ]
    ch_redux_dir_tumor  // channel: [mandatory] [ meta, redux_dir ]
    ch_tumor_rna_bam    // channel: [mandatory] [ meta, bam, bai ]

    // Reference data
    genome_fasta        // channel: [mandatory] /path/to/genome_fasta
    genome_version      // channel: [mandatory] genome version
    genome_fai          // channel: [mandatory] /path/to/genome_fai
    genome_dict         // channel: [mandatory] /path/to/genome_dict

    // Params
    sequencing_platform // string:  [mandatory] sequencing platform
    enable_germline     // boolean: [mandatory] Enable germline
    targeted_mode       // boolean: [mandatory] Set targeted mode

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    //
    // STEP: Handle inputs
    //
    def run_mode = Utils.getEnumFromString(params.mode, Constants.RunMode)
    def purity_estimate_mode = run_mode === Constants.RunMode.PURITY_ESTIMATE

    // Select input sources then sort
    // channel: runnable: [ meta, purple_dir, tumor_dna_bam, tumor_dna_bai, [redux_tsv_tumor, ...], tumor_rna_bam, tumor_rna_bai ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_purple_dir,
        ch_redux_dir_tumor,
        ch_tumor_rna_bam,
    )
        .map { meta, purple_dir, redux_dir_tumor, tumor_rna_bam, tumor_rna_bai ->

            def redux_dir_tumor_selected = Utils.selectCurrentOrExisting(redux_dir_tumor, meta, Constants.INPUT.REDUX_DIR_TUMOR)
            def (tumor_bam, tumor_bai) = Utils.getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)
            def redux_tsvs_tumor = Utils.getTumorReduxTsvs(meta, redux_dir_tumor_selected)

            return [
                meta,
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
                tumor_bam,
                tumor_bai,
                redux_tsvs_tumor,
                Utils.selectCurrentOrExisting(tumor_rna_bam, meta, Constants.INPUT.BAM_RNA_TUMOR),
                Utils.selectCurrentOrExisting(tumor_rna_bai, meta, Constants.INPUT.BAI_RNA_TUMOR),
            ]

        }
        .branch { meta, purple_dir, tumor_dna_bam, tumor_dna_bai, redux_tsvs_tumor, tumor_rna_bam, tumor_rna_bai ->
            def has_bam = tumor_dna_bam || tumor_rna_bam
            runnable: has_bam && purple_dir
            skip: true
                return meta
        }

    //
    // MODULE: SAGE append germline
    //
    // Select inputs that are eligible to run
    // channel: runnable: [ meta, purple_dir, tumor_dna_bam, tumor_dna_bai, [redux_tsv_tumor, ...], tumor_rna_bam, tumor_rna_bai ]
    // channel: skip: [ meta ]
    ch_inputs_germline_sorted = ch_inputs_sorted.runnable
        .branch { meta, purple_dir, tumor_dna_bam, tumor_dna_bai, redux_tsvs_tumor, tumor_rna_bam, tumor_rna_bai ->

            // NOTE(SW): explicit in expectation to always obtain the primary tumor DNA sample ID here
            def tumor_dna_id = Utils.getTumorDnaSampleName(meta, primary: true)

            def has_tumor_rna = Utils.hasTumorRna(meta)
            def has_normal_dna = Utils.hasNormalDna(meta)
            def has_smlv_germline = purple_dir.resolve("${tumor_dna_id}.purple.germline.vcf.gz").exists()

            def should_append_rna_variants = has_tumor_rna && has_normal_dna && has_smlv_germline

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SAGE_APPEND_DIR_NORMAL)

            runnable: should_append_rna_variants && ! has_existing && enable_germline
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_append, purple_smlv_vcf, [bam, ...], [bai, ...], [redux_tsvs_tumor, ...] ]
    ch_sage_append_germline_inputs = ch_inputs_germline_sorted.runnable
        .map { meta, purple_dir, tumor_dna_bam, tumor_dna_bai, redux_tsvs_tumor, tumor_rna_bam, tumor_rna_bai ->

            // NOTE(SW): explicit in expectation to always obtain the primary tumor DNA sample ID here
            def tumor_dna_id = Utils.getTumorDnaSampleName(meta, primary: true)
            def output_file_id = Utils.getNormalDnaSampleName(meta)

            def meta_append = [
                key: meta.group_id,
                id: "${meta.group_id}_${output_file_id}",
                output_file_id: output_file_id,
                reference_ids: [Utils.getTumorRnaSampleName(meta)],
            ]

            def bams = [tumor_rna_bam]
            def bais = [tumor_rna_bai]

            def purple_smlv_vcf = file(purple_dir).resolve("${tumor_dna_id}.purple.germline.vcf.gz")

            return [meta_append, purple_smlv_vcf, bams, bais, []]

        }

    // Run process
    SAGE_APPEND_GERMLINE(
        ch_sage_append_germline_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        sequencing_platform,
        targeted_mode,
    )

    ch_versions = ch_versions.mix(SAGE_APPEND_GERMLINE.out.versions)

    //
    // MODULE: SAGE append somatic
    //
    // Select inputs that are eligible to run
    // channel: runnable: [ meta, purple_dir, tumor_dna_bam, tumor_dna_bai, [redux_tsv_tumor, ...], tumor_rna_bam, tumor_rna_bai ]
    // channel: skip: [ meta ]
    ch_inputs_somatic_sorted = ch_inputs_sorted.runnable
        .branch { meta, purple_dir, tumor_dna_bam, tumor_dna_bai, redux_tsvs_tumor, tumor_rna_bam, tumor_rna_bai ->

            def tumor_dna_id = Utils.getTumorDnaSampleName(meta, primary: true)

            def has_tumor_rna = Utils.hasTumorRna(meta)
            def has_tumor_dna = Utils.hasTumorDna(meta)
            def has_smlv_somatic = file(purple_dir).resolve("${tumor_dna_id}.purple.somatic.vcf.gz")

            def should_append_rna_variants = ! purity_estimate_mode && has_tumor_rna && has_tumor_dna && has_smlv_somatic
            def should_append_longitudinal_variants = purity_estimate_mode && has_tumor_dna && has_smlv_somatic

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SAGE_APPEND_DIR_TUMOR)

            runnable: (should_append_rna_variants || should_append_longitudinal_variants) && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_append, purple_smlv_vcf, [bam, ...], [bai, ...], [redux_tsvs_tumor, ...] ]
    ch_sage_append_somatic_inputs = ch_inputs_somatic_sorted.runnable
        .map { meta, purple_dir, tumor_dna_bam, tumor_dna_bai, redux_tsvs_tumor, tumor_rna_bam, tumor_rna_bai ->

            // NOTE(SW): explicit in expectation to always obtain the primary tumor DNA sample ID here
            def tumor_dna_id = Utils.getTumorDnaSampleName(meta, primary: true)
            def output_file_id = purity_estimate_mode ? Utils.getTumorDnaSampleName(meta, primary: false) : tumor_dna_id

            def meta_append = [
                key: meta.group_id,
                id: "${meta.group_id}_${output_file_id}",
                output_file_id: output_file_id,
                reference_ids: [],
            ]

            def bams = []
            def bais = []
            def redux_tsvs = []

            if (! purity_estimate_mode && tumor_rna_bam) {
                meta_append.reference_ids.add(Utils.getTumorRnaSampleName(meta))
                bams.add(tumor_rna_bam)
                bais.add(tumor_rna_bai)
            }

            if (purity_estimate_mode && tumor_dna_bam) {
                meta_append.reference_ids.add(Utils.getTumorDnaSampleName(meta, primary: false))
                bams.add(tumor_dna_bam)
                bais.add(tumor_dna_bai)
                redux_tsvs = redux_tsvs_tumor
            }

            def purple_smlv_vcf = file(purple_dir).resolve("${tumor_dna_id}.purple.somatic.vcf.gz")

            return [meta_append, purple_smlv_vcf, bams, bais, redux_tsvs]
        }

    // Run process
    SAGE_APPEND_SOMATIC(
        ch_sage_append_somatic_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        sequencing_platform,
        targeted_mode,
    )

    ch_versions = ch_versions.mix(SAGE_APPEND_SOMATIC.out.versions)

    //
    // STEP: Handle outputs
    //
    // Set outputs, restoring original meta
    // channel: [ meta, sage_append_dir ]
    ch_outputs_somatic = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(SAGE_APPEND_SOMATIC.out.sage_append_dir, ch_inputs),
            ch_inputs_somatic_sorted.skip.map { meta -> [meta, []] },
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    // channel: [ meta, sage_append_dir ]
    ch_outputs_germline = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(SAGE_APPEND_GERMLINE.out.sage_append_dir, ch_inputs),
            ch_inputs_germline_sorted.skip.map { meta -> [meta, []] },
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    somatic_dir  = ch_outputs_somatic  // channel: [ meta, sage_append_dir ]
    germline_dir = ch_outputs_germline // channel: [ meta, sage_append_dir ]

    versions     = ch_versions         // channel: [ versions.yml ]
}
