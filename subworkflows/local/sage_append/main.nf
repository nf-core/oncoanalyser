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
    ch_inputs         // channel: [mandatory] [ meta ]
    ch_purple_dir     // channel: [mandatory] [ meta, purple_dir ]
    ch_tumor_dna_bam  // channel: [mandatory] [ meta, bam, bai ]
    ch_tumor_dna_tsv  // channel: [mandatory] [ meta, dup_freq_tsv, jitter_tsv, ms_tsv ]
    ch_tumor_rna_bam  // channel: [mandatory] [ meta, bam, bai ]

    // Reference data
    genome_fasta     // channel: [mandatory] /path/to/genome_fasta
    genome_version   // channel: [mandatory] genome version
    genome_fai       // channel: [mandatory] /path/to/genome_fai
    genome_dict      // channel: [mandatory] /path/to/genome_dict

    // Params
    enable_germline  // boolean: [mandatory] Enable germline
    is_targeted_mode // boolean: [mandatory] Running in targeted/panel mode?

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources and sort
    // channel: runnable: [ meta, tumor_dna_bam, tumor_dna_bai, [tumor_dna_redux_tsv, ...], tumor_rna_bam, tumor_rna_bai, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_tumor_dna_bam,
        ch_tumor_dna_tsv,
        ch_tumor_rna_bam,
        ch_purple_dir,
    )
        .map { meta, tumor_dna_bam, tumor_dna_bai, tumor_dna_dup_freq_tsv, tumor_dna_jitter_tsv, tumor_dna_ms_tsv, tumor_rna_bam, tumor_rna_bai, purple_dir ->

            def tumor_dna_redux_tsv_list = [
                tumor_dna_jitter_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_JITTER_TSV_TUMOR),
                tumor_dna_ms_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_MS_TSV_TUMOR),
            ]

            tumor_dna_redux_tsv_list = tumor_dna_redux_tsv_list.findAll{ it != [] }

            return [
                meta,
                Utils.selectCurrentOrExisting(tumor_dna_bam, meta, Constants.INPUT.BAM_REDUX_DNA_TUMOR),
                tumor_dna_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_TUMOR),
                tumor_dna_redux_tsv_list,
                Utils.selectCurrentOrExisting(tumor_rna_bam, meta, Constants.INPUT.BAM_RNA_TUMOR),
                tumor_rna_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_RNA_TUMOR),
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
            ]
        }
        .branch { meta, tumor_dna_bam, tumor_dna_bai, tumor_dna_redux_tsv, tumor_rna_bam, tumor_rna_bai, purple_dir ->
            def has_bam = tumor_dna_bam || tumor_rna_bam
            runnable: has_bam && purple_dir
            skip: true
                return meta
        }

    //
    // MODULE: SAGE append germline
    //
    // Select inputs that are eligible to run
    // channel: runnable: [ meta, tumor_dna_bam, tumor_dna_bai, [tumor_dna_redux_tsv, ...], tumor_rna_bam, tumor_rna_bai, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_germline_sorted = ch_inputs_sorted.runnable
        .branch { meta, tumor_dna_bam, tumor_dna_bai, tumor_dna_redux_tsv, tumor_rna_bam, tumor_rna_bai, purple_dir ->

            // NOTE(SW): explicit in expectation to always obtain the primary tumor DNA sample ID here
            def tumor_dna_id = Utils.getTumorDnaSampleName(meta, primary: true)

            def has_normal_dna = Utils.hasNormalDna(meta)
            def has_tumor_rna = Utils.hasTumorRna(meta)

            def has_smlv_germline = file(purple_dir).resolve("${tumor_dna_id}.purple.germline.vcf.gz")
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SAGE_APPEND_DIR_NORMAL)

            runnable: has_normal_dna && has_tumor_rna && has_smlv_germline && !has_existing && enable_germline
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_append, purple_smlv_vcf, [bam, ...], [bai, ...], [tumor_dna_redux_tsv, ...] ]
    ch_sage_append_germline_inputs = ch_inputs_germline_sorted.runnable
        .map { meta, tumor_dna_bam, tumor_dna_bai, tumor_dna_redux_tsv, tumor_rna_bam, tumor_rna_bai, purple_dir ->

            // NOTE(SW): explicit in expectation to always obtain the primary tumor DNA sample ID here
            def tumor_dna_id   = Utils.getTumorDnaSampleName(meta, primary: true)
            def output_file_id = Utils.getNormalDnaSampleName(meta)

            def meta_append = [
                key: meta.group_id,
                id: meta.group_id,
                output_file_id: output_file_id,
                reference_ids: [Utils.getTumorRnaSampleName(meta)],
            ]

            def bams = [tumor_rna_bam]
            def bais = [tumor_rna_bai]

            def purple_smlv_vcf = file(purple_dir).resolve("${tumor_dna_id}.purple.germline.vcf.gz")

            return [meta_append, purple_smlv_vcf, bams, bais, tumor_dna_redux_tsv]
        }

    // Run process
    SAGE_APPEND_GERMLINE(
        ch_sage_append_germline_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        is_targeted_mode,
    )

    ch_versions = ch_versions.mix(SAGE_APPEND_GERMLINE.out.versions)

    //
    // MODULE: SAGE append somatic
    //
    // Select inputs that are eligible to run
    // channel: runnable: [ meta, tumor_dna_bam, tumor_dna_bai, [tumor_dna_redux_tsv, ...], tumor_rna_bam, tumor_rna_bai, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_somatic_sorted = ch_inputs_sorted.runnable
        .branch { meta, tumor_dna_bam, tumor_dna_bai, tumor_dna_redux_tsv, tumor_rna_bam, tumor_rna_bai, purple_dir ->

            def tumor_dna_id = Utils.getTumorDnaSampleName(meta, primary: true)

            def has_tumor_sequence_data = Utils.hasTumorRna(meta) || Utils.hasTumorDna(meta)

            def has_smlv_somatic = file(purple_dir).resolve("${tumor_dna_id}.purple.somatic.vcf.gz")
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SAGE_APPEND_DIR_TUMOR)

            runnable: has_tumor_sequence_data && has_smlv_somatic && !has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_append, purple_smlv_vcf, [bam, ...], [bai, ...], [tumor_dna_redux_tsv, ...] ]
    ch_sage_append_somatic_inputs = ch_inputs_somatic_sorted.runnable
        .map { meta, tumor_dna_bam, tumor_dna_bai, tumor_dna_redux_tsv, tumor_rna_bam, tumor_rna_bai, purple_dir ->

            def run_mode = Utils.getEnumFromString(params.mode, Constants.RunMode)

            def tumor_dna_id = Utils.getTumorDnaSampleName(meta, primary: true)

            def output_file_id = run_mode === Constants.RunMode.PURITY_ESTIMATE
                ? Utils.getTumorDnaSampleName(meta, primary: false)
                : tumor_dna_id

            def meta_append = [
                key: meta.group_id,
                id: meta.group_id,
                output_file_id: output_file_id,
                reference_ids: [],
            ]

            def bams = []
            def bais = []

            if (tumor_dna_bam) {
                meta_append.reference_ids.add(Utils.getTumorDnaSampleName(meta))
                bams.add(tumor_dna_bam)
                bais.add(tumor_dna_bai)
            }

            if (tumor_rna_bam) {
                meta_append.reference_ids.add(Utils.getTumorRnaSampleName(meta))
                bams.add(tumor_rna_bam)
                bais.add(tumor_rna_bai)
            }

            def purple_smlv_vcf = file(purple_dir).resolve("${tumor_dna_id}.purple.somatic.vcf.gz")

            return [meta_append, purple_smlv_vcf, bams, bais, tumor_dna_redux_tsv]
        }

    // Run process
    SAGE_APPEND_SOMATIC(
        ch_sage_append_somatic_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        is_targeted_mode,
    )

    ch_versions = ch_versions.mix(SAGE_APPEND_SOMATIC.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, sage_append_dir ]
    ch_somatic_dir = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(SAGE_APPEND_SOMATIC.out.sage_append_dir, ch_inputs),
            ch_inputs_somatic_sorted.skip.map { meta -> [meta, []] },
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    ch_germline_dir = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(SAGE_APPEND_GERMLINE.out.sage_append_dir, ch_inputs),
            ch_inputs_germline_sorted.skip.map { meta -> [meta, []] },
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    somatic_dir  = ch_somatic_dir  // channel: [ meta, sage_append_dir ]
    germline_dir = ch_germline_dir // channel: [ meta, sage_append_dir ]

    versions     = ch_versions     // channel: [ versions.yml ]
}
