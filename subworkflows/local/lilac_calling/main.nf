//
// LILAC is a WGS tool for HLA typing and somatic CNV and SNV calling
//

include { LILAC } from '../../../modules/local/lilac/main'

workflow LILAC_CALLING {
    take:
    // Sample data
    ch_inputs           // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor  // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal // channel: [mandatory] [ meta, redux_dir ]
    ch_tumor_rna_aln    // channel: [mandatory] [ meta, aln, idx ]
    ch_purple           // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_fasta        // channel: [mandatory] /path/to/genome_fasta
    genome_version      // channel: [mandatory] genome version
    genome_fai          // channel: [mandatory] /path/to/genome_fai
    lilac_resource_dir  // channel: [mandatory] /path/to/lilac_resource_dir/

    // Params
    sequencing_platform // string:  [mandatory] sequencing platform
    targeted_mode       // boolean: [mandatory] Set targeted mode

    main:
    // Select input sources then sort
    // channel: runnable: [meta, normal_dna_aln, normal_dna_idx, tumor_dna_aln, tumor_dna_idx, tumor_rna_aln, tumor_rna_idx, purple_dir]
    // channel: skip: [ meta ]
    ch_dna_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_redux_dir_tumor,
        ch_redux_dir_normal,
        ch_tumor_rna_aln,
        ch_purple,
    )
        .map { meta, redux_dir_tumor, redux_dir_normal, tumor_rna_aln, tumor_rna_idx, purple_dir ->

            def redux_dir_tumor_selected = Utils.selectCurrentOrExisting(redux_dir_tumor, meta, Constants.INPUT.REDUX_DIR_TUMOR)
            def redux_dir_normal_selected = Utils.selectCurrentOrExisting(redux_dir_normal, meta, Constants.INPUT.REDUX_DIR_NORMAL)

            def (tumor_dna_aln, tumor_dna_idx) = Utils.getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)
            def (normal_dna_aln, normal_dna_idx) = Utils.getNormalReduxDirAlignment(meta, redux_dir_normal_selected)

            return [
                meta,
                normal_dna_aln,
                normal_dna_idx,
                tumor_dna_aln,
                tumor_dna_idx,
                Utils.selectCurrentOrExisting(tumor_rna_aln, meta, Constants.INPUT.ALN_RNA_TUMOR),
                Utils.selectCurrentOrExisting(tumor_rna_idx, meta, Constants.INPUT.IDX_RNA_TUMOR),
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
            ]

        }
        .branch { meta, normal_dna_aln, normal_dna_idx, tumor_dna_aln, tumor_dna_idx, tumor_rna_aln, tumor_rna_idx, purple_dir ->

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.LILAC_DIR)

            def tumor_normal_mode = tumor_dna_aln && normal_dna_aln

            def tumor_dna_id = Utils.getTumorDnaSampleName(meta)
            def has_tn_smlv_vcf = purple_dir ? purple_dir.resolve("${tumor_dna_id}.purple.somatic.vcf.gz").exists() : false

            runnable: (tumor_dna_aln || normal_dna_aln) && (has_tn_smlv_vcf || ! tumor_normal_mode) && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_lilac, normal_dna_aln, normal_dna_idx, tumor_dna_aln, tumor_dna_idx, tumor_rna_aln, tumor_rna_idx, purple_dir
    ch_lilac_inputs = WorkflowOncoanalyser.groupByMeta(
        ch_dna_inputs_sorted.runnable,
    )
        .map { meta, normal_dna_aln, normal_dna_idx, tumor_dna_aln, tumor_dna_idx, tumor_rna_aln, tumor_rna_idx, purple_dir ->

            def meta_lilac = [
                key: meta.group_id,
                id: meta.group_id,
            ]

            if (Utils.hasTumorDna(meta)) {
                meta_lilac.tumor_id = Utils.getTumorDnaSampleName(meta)
            }

            if (Utils.hasNormalDna(meta)) {
                meta_lilac.normal_id = Utils.getNormalDnaSampleName(meta)
            }

            return [meta_lilac, normal_dna_aln, normal_dna_idx, tumor_dna_aln, tumor_dna_idx, tumor_rna_aln, tumor_rna_idx, purple_dir]
        }

    // Run process
    LILAC(
        ch_lilac_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        lilac_resource_dir,
        targeted_mode,
        sequencing_platform,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, lilac_dir ]
    ch_outputs = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('lilac_dir'), ch_inputs),
            ch_dna_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    lilac_dir = ch_outputs // channel: [ meta, lilac_dir ]
}
