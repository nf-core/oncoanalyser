//
// SAGE append adds additional sample data to an existing SAGE VCF
//

include { SAGE_APPEND as SAGE_APPEND_SOMATIC  } from '../../../modules/local/sage/append/main'
include { SAGE_APPEND as SAGE_APPEND_GERMLINE } from '../../../modules/local/sage/append/main'

workflow SAGE_APPEND {
    take:
    // Sample data
    ch_inputs           // channel: [mandatory] [ meta ]
    ch_purple_dir       // channel: [mandatory] [ meta, purple_dir ]
    ch_redux_dir_tumor  // channel: [mandatory] [ meta, redux_dir ]
    ch_tumor_rna_aln    // channel: [mandatory] [ meta, aln, idx ]

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
    //
    // STEP: Handle inputs
    //
    def run_mode = Utils.getEnumFromString(params.mode, Constants.RunMode)
    def purity_estimate_mode = run_mode == Constants.RunMode.PURITY_ESTIMATE

    // Select input sources then sort
    // channel: runnable: [ meta, purple_dir, tumor_dna_aln, tumor_dna_idx, [redux_tsv_tumor, ...], tumor_rna_aln, tumor_rna_idx ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_purple_dir,
        ch_redux_dir_tumor,
        ch_tumor_rna_aln,
    )
        .map { meta, purple_dir, redux_dir_tumor, tumor_rna_aln, tumor_rna_idx ->

            def redux_dir_tumor_selected = Utils.selectCurrentOrExisting(redux_dir_tumor, meta, Constants.INPUT.REDUX_DIR_TUMOR)
            def (tumor_aln, tumor_idx) = Utils.getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)
            def redux_tsvs_tumor = Utils.getTumorReduxTsvs(meta, redux_dir_tumor_selected)

            return [
                meta,
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
                tumor_aln,
                tumor_idx,
                redux_tsvs_tumor,
                Utils.selectCurrentOrExisting(tumor_rna_aln, meta, Constants.INPUT.ALN_RNA_TUMOR),
                Utils.selectCurrentOrExisting(tumor_rna_idx, meta, Constants.INPUT.IDX_RNA_TUMOR),
            ]

        }
        .branch { meta, purple_dir, tumor_dna_aln, tumor_dna_idx, redux_tsvs_tumor, tumor_rna_aln, tumor_rna_idx ->
            def has_aln = tumor_dna_aln || tumor_rna_aln
            runnable: has_aln && purple_dir
            skip: true
                return meta
        }

    //
    // MODULE: SAGE append germline
    //
    // Select inputs that are eligible to run
    // channel: runnable: [ meta, purple_dir, tumor_dna_aln, tumor_dna_idx, [redux_tsv_tumor, ...], tumor_rna_aln, tumor_rna_idx ]
    // channel: skip: [ meta ]
    ch_inputs_germline_sorted = ch_inputs_sorted.runnable
        .branch { meta, purple_dir, tumor_dna_aln, tumor_dna_idx, redux_tsvs_tumor, tumor_rna_aln, tumor_rna_idx ->

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
    // channel: [ meta_append, purple_smlv_vcf, [aln, ...], [idx, ...], [redux_tsvs_tumor, ...] ]
    ch_sage_append_germline_inputs = ch_inputs_germline_sorted.runnable
        .map { meta, purple_dir, _tumor_dna_aln, _tumor_dna_idx, _redux_tsvs_tumor, tumor_rna_aln, tumor_rna_idx ->

            // NOTE(SW): explicit in expectation to always obtain the primary tumor DNA sample ID here
            def tumor_dna_id = Utils.getTumorDnaSampleName(meta, primary: true)
            def output_file_id = Utils.getNormalDnaSampleName(meta)

            def meta_append = [
                key: meta.group_id,
                topic_key: 'germline',
                id: "${meta.group_id}_${output_file_id}",
                output_file_id: output_file_id,
                reference_ids: [Utils.getTumorRnaSampleName(meta)],
            ]

            def alns = [tumor_rna_aln]
            def idxs = [tumor_rna_idx]

            def purple_smlv_vcf = purple_dir.resolve("${tumor_dna_id}.purple.germline.vcf.gz")

            // NOTE(SW): in stub mode we allow absence of indexes so must handle here
            idxs = idxs.flatten()

            return [meta_append, purple_smlv_vcf, alns, idxs, []]

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

    //
    // MODULE: SAGE append somatic
    //
    // Select inputs that are eligible to run
    // channel: runnable: [ meta, purple_dir, tumor_dna_aln, tumor_dna_idx, [redux_tsv_tumor, ...], tumor_rna_aln, tumor_rna_idx ]
    // channel: skip: [ meta ]
    ch_inputs_somatic_sorted = ch_inputs_sorted.runnable
        .branch { meta, purple_dir, tumor_dna_aln, tumor_dna_idx, redux_tsvs_tumor, tumor_rna_aln, tumor_rna_idx ->

            def tumor_dna_id = Utils.getTumorDnaSampleName(meta, primary: true)

            def has_tumor_rna = Utils.hasTumorRna(meta)
            def has_tumor_dna = Utils.hasTumorDna(meta)
            def has_smlv_somatic = purple_dir.resolve("${tumor_dna_id}.purple.somatic.vcf.gz").exists()

            def should_append_rna_variants = ! purity_estimate_mode && has_tumor_rna && has_tumor_dna && has_smlv_somatic
            def should_append_longitudinal_variants = purity_estimate_mode && has_tumor_dna && has_smlv_somatic

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SAGE_APPEND_DIR_TUMOR)

            runnable: (should_append_rna_variants || should_append_longitudinal_variants) && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_append, purple_smlv_vcf, [aln, ...], [idx, ...], [redux_tsvs_tumor, ...] ]
    ch_sage_append_somatic_inputs = ch_inputs_somatic_sorted.runnable
        .map { meta, purple_dir, tumor_dna_aln, tumor_dna_idx, redux_tsvs_tumor, tumor_rna_aln, tumor_rna_idx ->

            // NOTE(SW): explicit in expectation to always obtain the primary tumor DNA sample ID here
            def tumor_dna_id = Utils.getTumorDnaSampleName(meta, primary: true)
            def output_file_id = purity_estimate_mode ? Utils.getTumorDnaSampleName(meta, primary: false) : tumor_dna_id

            def meta_append = [
                key: meta.group_id,
                topic_key: 'somatic',
                id: "${meta.group_id}_${output_file_id}",
                output_file_id: output_file_id,
                reference_ids: [],
            ]

            def alns = []
            def idxs = []
            def redux_tsvs = []

            if (! purity_estimate_mode && tumor_rna_aln) {
                meta_append.reference_ids.add(Utils.getTumorRnaSampleName(meta))
                alns.add(tumor_rna_aln)
                idxs.add(tumor_rna_idx)
            }

            if (purity_estimate_mode && tumor_dna_aln) {
                meta_append.reference_ids.add(Utils.getTumorDnaSampleName(meta, primary: false))
                alns.add(tumor_dna_aln)
                idxs.add(tumor_dna_idx)
                redux_tsvs = redux_tsvs_tumor
            }

            def purple_smlv_vcf = purple_dir.resolve("${tumor_dna_id}.purple.somatic.vcf.gz")

            // NOTE(SW): in stub mode we allow absence of indexes so must handle here
            idxs = idxs.flatten()

            return [meta_append, purple_smlv_vcf, alns, idxs, redux_tsvs]
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

    //
    // STEP: Handle outputs
    //
    // Set outputs, restoring original meta
    // channel: [ meta, sage_append_dir ]
    ch_outputs_somatic = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('sage_append_dir').filter { d -> d[0].topic_key == 'somatic' }, ch_inputs),
            ch_inputs_somatic_sorted.skip.map { meta -> [meta, []] },
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    // channel: [ meta, sage_append_dir ]
    ch_outputs_germline = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('sage_append_dir').filter { d -> d[0].topic_key == 'germline' }, ch_inputs),
            ch_inputs_germline_sorted.skip.map { meta -> [meta, []] },
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    somatic_dir  = ch_outputs_somatic  // channel: [ meta, sage_append_dir ]
    germline_dir = ch_outputs_germline // channel: [ meta, sage_append_dir ]
}
