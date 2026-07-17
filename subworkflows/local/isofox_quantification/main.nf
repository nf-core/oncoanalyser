//
// Isofox estimates transcript abundance, detects novel SJs, and identifies fusion events
//

include { ISOFOX } from '../../../modules/local/isofox/run/main'

workflow ISOFOX_QUANTIFICATION {
    take:
    // Sample data
    ch_inputs                  // channel: [mandatory] [ meta ]
    ch_tumor_rna_aln           // channel: [mandatory] [ meta, aln, idx ]

    // Reference data
    genome_fasta               // channel: [mandatory] /path/to/genome_fasta
    genome_version             // channel: [mandatory] genome version
    genome_fai                 // channel: [mandatory] /path/to/genome_fai
    ensembl_data_resources     // channel: [mandatory] /path/to/ensembl_data_resources/
    driver_gene_panel          // channel: [mandatory] /path/to/driver_gene_panel
    known_fusion_data          // channel: [mandatory] /path/to/known_fusion_data
    isofox_excluded_regions    // channel: [mandatory] /path/to/isofox_excluded_regions
    isofox_gene_distribution   // channel: [mandatory] /path/to/isofox_gene_distribution
    isofox_alt_sj_distribution // channel: [mandatory] /path/to/isofox_alt_sj_distribution
    isofox_counts              // channel: [mandatory] /path/to/isofox_counts
    isofox_gc_ratios           // channel: [mandatory] /path/to/isofox_gc_ratios
    isofox_tpm_norm            // channel: [optional]  /path/to/isofox_tpm_norm

    // Params
    isofox_functions       //  string: [optional]  Isofox functions
    isofox_read_length     //  string: [mandatory] Isofox read length

    main:
    // Select input sources then sort
    // channel: runnable: [ meta, tumor_aln, tumor_idx ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_tumor_rna_aln
        .map { meta, tumor_aln, tumor_idx ->
            return [
                meta,
                Utils.selectCurrentOrExisting(tumor_aln, meta, Constants.INPUT.ALN_RNA_TUMOR),
                Utils.selectCurrentOrExisting(tumor_idx, meta, Constants.INPUT.IDX_RNA_TUMOR),
            ]
        }
        .branch { meta, tumor_aln, tumor_idx ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.ISOFOX_DIR)
            runnable: tumor_aln && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_isofox, tumor_aln, tumor_idx ]
    ch_isofox_inputs = ch_inputs_sorted.runnable
        .map { meta, tumor_aln, tumor_idx ->

            def meta_isofox = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta) ?: Utils.getTumorRnaSampleName(meta),
            ]

            return [meta_isofox, tumor_aln, tumor_idx]
        }

    // Run process
    ISOFOX(
        ch_isofox_inputs,
        isofox_functions,
        isofox_read_length,
        genome_fasta,
        genome_version,
        genome_fai,
        ensembl_data_resources,
        driver_gene_panel,
        known_fusion_data,
        isofox_excluded_regions,
        isofox_gene_distribution,
        isofox_alt_sj_distribution,
        isofox_counts,
        isofox_gc_ratios,
        isofox_tpm_norm,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, isofox_dir ]
    ch_outputs = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('isofox_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    isofox_dir = ch_outputs // channel: [ meta, isofox_dir ]
}
