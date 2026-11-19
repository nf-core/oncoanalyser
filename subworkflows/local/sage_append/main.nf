//
// SAGE append adds data from an additional tumor sample to an existing Purple somatic VCF
//
// The Purple somatic VCF always belongs to the primary tumor DNA sample. The appended sample is selected by the caller:
//   - wgts/targeted: tumor RNA sample
//   - purity estimate: longitudinal tumor DNA sample
//

include { SAGE_APPEND as SAGE_APPEND_SOMATIC } from '../../../modules/local/sage/append/main'

workflow SAGE_APPEND {
    take:
    // Sample data
    ch_inputs            // channel: [mandatory] [ meta ]
    ch_purple_dir        // channel: [mandatory] [ meta, purple_dir ] (unused in longitudinal mode)
    ch_redux_dir         // channel: [mandatory] [ meta, redux_dir ] for the appended sample

    // Reference data
    genome_fasta         // channel: [mandatory] /path/to/genome_fasta
    genome_version       // channel: [mandatory] genome version
    genome_fai           // channel: [mandatory] /path/to/genome_fai
    genome_dict          // channel: [mandatory] /path/to/genome_dict

    // Params
    sequencing_platform  // string:  [mandatory] sequencing platform
    targeted_mode        // boolean: [mandatory] Set targeted mode
    append_sample_type   // string:  [mandatory] 'rna' (tumor RNA) or 'longitudinal' (longitudinal tumor DNA)

    main:

    // Select input sources then sort
    // channel: runnable: [ meta, append_sample_id, output_file_id, purple_smlv_vcf, aln, idx, [redux_tsv, ...] ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_purple_dir,
        ch_redux_dir,
    )
        .map { meta, purple_dir, redux_dir ->

            def primary_tumor_dna_id = Utils.getTumorDnaSampleName(meta, primary: true)

            def purple_dir_selected
            def append_sample_id
            def output_file_id
            def redux_dir_selected
            def redux_tsvs

            if (append_sample_type == 'rna') {

                // Purple dir from this pipeline run or samplesheet
                purple_dir_selected = Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR)

                def tumor_rna_id = Utils.getTumorRnaSampleName(meta)
                redux_dir_selected = Utils.selectCurrentOrExisting(redux_dir, meta, Constants.INPUT.REDUX_DIR_RNA)
                redux_tsvs = []

                append_sample_id = tumor_rna_id
                output_file_id   = primary_tumor_dna_id

            } else if (append_sample_type == 'longitudinal') {

                // NOTE(LN): The primary sample Purple dir is always from the samplesheet (i.e. from a previous pipeline run)
                purple_dir_selected = Utils.getInput(meta, Constants.INPUT.PURPLE_DIR)

                def longitudinal_tumor_dna_id = Utils.getTumorDnaSampleName(meta, primary: false)
                redux_dir_selected = Utils.selectCurrentOrExisting(redux_dir, meta, Constants.INPUT.REDUX_DIR_TUMOR)
                redux_tsvs = Utils.getReduxTsvs(longitudinal_tumor_dna_id, redux_dir_selected)

                append_sample_id = longitudinal_tumor_dna_id
                output_file_id   = longitudinal_tumor_dna_id

            } else {
                error("Invalid Sage append sample type '${append_sample_type}', expected 'rna' or 'longitudinal'")
            }

            def purple_smlv_vcf = purple_dir_selected
                ? purple_dir_selected.resolve("${primary_tumor_dna_id}.purple.somatic.vcf.gz")
                : []

            def (aln, idx) = Utils.getReduxDirAlignment(append_sample_id, redux_dir_selected)

            return [meta, append_sample_id, output_file_id, purple_smlv_vcf, aln, idx, redux_tsvs]
        }
        .branch { meta, append_sample_id, output_file_id, purple_smlv_vcf, aln, idx, redux_tsvs ->

            def has_tumor_dna = Utils.hasTumorDna(meta)
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SAGE_APPEND_DIR_TUMOR)

            runnable: has_tumor_dna && aln && purple_smlv_vcf && purple_smlv_vcf.exists() && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_append, purple_smlv_vcf, [aln, ...], [idx, ...], [redux_tsv, ...] ]
    ch_sage_append_inputs = ch_inputs_sorted.runnable
        .map { meta, append_sample_id, output_file_id, purple_smlv_vcf, aln, idx, redux_tsvs ->

            def meta_append = [
                key: meta.group_id,
                id: "${meta.group_id}_${output_file_id}",
                output_file_id: output_file_id,
                reference_ids: [append_sample_id],
            ]

            // NOTE(SW): in stub mode we allow absence of indexes so must handle here
            def idxs = [idx].flatten()

            return [meta_append, purple_smlv_vcf, [aln], idxs, redux_tsvs]
        }

    // Run process
    SAGE_APPEND_SOMATIC(
        ch_sage_append_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        sequencing_platform,
        targeted_mode,
        append_sample_type == 'longitudinal',  // write_frag_lengths
    )

    // Set outputs, restoring original meta
    // channel: [ meta, sage_append_dir ]
    ch_outputs = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('sage_append_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    sage_append_dir = ch_outputs // channel: [ meta, sage_append_dir ]
}
