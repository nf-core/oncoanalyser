//
// SAGE append adds longitudinal tumor DNA sample data to an existing Purple somatic VCF
//

include { SAGE_APPEND } from '../../../modules/local/sage/append/main'

workflow SAGE_APPEND_LONGITUDINAL {
    take:
    // Sample data
    ch_inputs           // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor  // channel: [mandatory] [ meta, redux_dir ]

    // Reference data
    genome_fasta        // channel: [mandatory] /path/to/genome_fasta
    genome_version      // channel: [mandatory] genome version
    genome_fai          // channel: [mandatory] /path/to/genome_fai
    genome_dict         // channel: [mandatory] /path/to/genome_dict

    // Params
    sequencing_platform // string:  [mandatory] sequencing platform
    targeted_mode       // boolean: [mandatory] Set targeted mode

    main:
    // Select input sources then sort
    // channel: runnable: [ meta, purple_dir, tumor_dna_aln, tumor_dna_idx, [redux_tsv_tumor, ...] ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_redux_dir_tumor
        .map { meta, redux_dir_tumor ->

            // NOTE(LN): PURPLE is not run in purity estimate mode, the primary sample PURPLE directory is always an existing input
            def purple_dir = Utils.getInput(meta, Constants.INPUT.PURPLE_DIR)

            def redux_dir_tumor_selected = Utils.selectCurrentOrExisting(redux_dir_tumor, meta, Constants.INPUT.REDUX_DIR_TUMOR)
            def (tumor_dna_aln, tumor_dna_idx) = Utils.getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)
            def redux_tsvs_tumor = Utils.getTumorReduxTsvs(meta, redux_dir_tumor_selected)

            return [
                meta,
                purple_dir,
                tumor_dna_aln,
                tumor_dna_idx,
                redux_tsvs_tumor,
            ]

        }
        .branch { meta, purple_dir, tumor_dna_aln, tumor_dna_idx, redux_tsvs_tumor ->

            def tumor_dna_id = Utils.getTumorDnaSampleName(meta, primary: true)

            def has_tumor_dna = Utils.hasTumorDna(meta)
            def has_smlv_somatic = purple_dir ? purple_dir.resolve("${tumor_dna_id}.purple.somatic.vcf.gz").exists() : false

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SAGE_APPEND_DIR_TUMOR)

            runnable: tumor_dna_aln && has_tumor_dna && has_smlv_somatic && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_append, purple_smlv_vcf, [aln, ...], [idx, ...], [redux_tsvs_tumor, ...] ]
    ch_sage_append_inputs = ch_inputs_sorted.runnable
        .map { meta, purple_dir, tumor_dna_aln, tumor_dna_idx, redux_tsvs_tumor ->

            def tumor_dna_id = Utils.getTumorDnaSampleName(meta, primary: true)
            def longitudinal_id = Utils.getTumorDnaSampleName(meta, primary: false)

            def meta_append = [
                key: meta.group_id,
                topic_key: 'longitudinal',
                id: "${meta.group_id}_${longitudinal_id}",
                output_file_id: longitudinal_id,
                reference_ids: [longitudinal_id],
            ]

            def purple_smlv_vcf = purple_dir.resolve("${tumor_dna_id}.purple.somatic.vcf.gz")

            // NOTE(SW): in stub mode we allow absence of indexes so must handle here
            def idxs = [tumor_dna_idx].flatten()

            return [meta_append, purple_smlv_vcf, [tumor_dna_aln], idxs, redux_tsvs_tumor]
        }

    // Run process
    SAGE_APPEND(
        ch_sage_append_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        sequencing_platform,
        targeted_mode,
        true,  // write_frag_lengths
    )

    // Set outputs, restoring original meta
    // channel: [ meta, sage_append_dir ]
    ch_outputs = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('sage_append_dir').filter { d -> d[0].topic_key == 'longitudinal' }, ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    sage_append_dir = ch_outputs // channel: [ meta, sage_append_dir ]
}
