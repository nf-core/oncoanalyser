//
// SAGE append adds tumor RNA sample data to an existing Purple somatic VCF
//

include { SAGE_APPEND } from '../../../modules/local/sage/append/main'

workflow SAGE_APPEND_RNA {
    take:
    // Sample data
    ch_inputs           // channel: [mandatory] [ meta ]
    ch_purple_dir       // channel: [mandatory] [ meta, purple_dir ]
    ch_redux_dir_rna    // channel: [mandatory] [ meta, redux_dir ]

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
    // channel: runnable: [ meta, purple_dir, tumor_rna_aln, tumor_rna_idx ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_purple_dir,
        ch_redux_dir_rna,
    )
        .map { meta, purple_dir, redux_dir_rna ->

            def purple_dir_selected = Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR)
            def redux_dir_rna_selected = Utils.selectCurrentOrExisting(redux_dir_rna, meta, Constants.INPUT.REDUX_DIR_RNA)
            def (tumor_rna_aln, tumor_rna_idx) = Utils.getTumorRnaReduxDirAlignment(meta, redux_dir_rna_selected)

            return [
                meta,
                purple_dir_selected,
                tumor_rna_aln,
                tumor_rna_idx,
            ]

        }
        .branch { meta, purple_dir, tumor_rna_aln, tumor_rna_idx ->

            def tumor_dna_id = Utils.getTumorDnaSampleName(meta, primary: true)

            def has_tumor_rna = Utils.hasTumorRna(meta)
            def has_tumor_dna = Utils.hasTumorDna(meta)
            def has_smlv_somatic = purple_dir ? purple_dir.resolve("${tumor_dna_id}.purple.somatic.vcf.gz").exists() : false

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SAGE_APPEND_DIR_TUMOR)

            runnable: tumor_rna_aln && has_tumor_rna && has_tumor_dna && has_smlv_somatic && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_append, purple_smlv_vcf, [aln, ...], [idx, ...], [redux_tsvs_tumor, ...] ]
    ch_sage_append_inputs = ch_inputs_sorted.runnable
        .map { meta, purple_dir, tumor_rna_aln, tumor_rna_idx ->

            def tumor_dna_id = Utils.getTumorDnaSampleName(meta, primary: true)

            def meta_append = [
                key: meta.group_id,
                topic_key: 'rna',
                id: "${meta.group_id}_${tumor_dna_id}",
                output_file_id: tumor_dna_id,
                reference_ids: [Utils.getTumorRnaSampleName(meta)],
            ]

            def purple_smlv_vcf = purple_dir.resolve("${tumor_dna_id}.purple.somatic.vcf.gz")

            // NOTE(SW): in stub mode we allow absence of indexes so must handle here
            def idxs = [tumor_rna_idx].flatten()

            return [meta_append, purple_smlv_vcf, [tumor_rna_aln], idxs, []]
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
        false,  // write_frag_lengths
    )

    // Set outputs, restoring original meta
    // channel: [ meta, sage_append_dir ]
    ch_outputs = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('sage_append_dir').filter { d -> d[0].topic_key == 'rna' }, ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    sage_append_dir = ch_outputs // channel: [ meta, sage_append_dir ]
}
