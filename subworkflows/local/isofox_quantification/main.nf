//
// Isofox estimates transcript abundance, detects novel SJs, and identifies fusion events
//

import Constants
import Utils

include { ISOFOX } from '../../../modules/local/isofox/main'

workflow ISOFOX_QUANTIFICATION {
    take:
    // Sample data
    ch_inputs              // channel: [mandatory] [ meta ]
    ch_tumor_rna_bam       // channel: [mandatory] [ meta, bam, bai ]

    // Reference data
    genome_fasta           // channel: [mandatory] /path/to/genome_fasta
    genome_version         // channel: [mandatory] genome version
    genome_fai             // channel: [mandatory] /path/to/genome_fai
    ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/
    known_fusion_data      // channel: [mandatory] /path/to/known_fusion_data
    isofox_counts          // channel: [mandatory] /path/to/isofox_counts
    isofox_gc_ratios       // channel: [mandatory] /path/to/isofox_gc_ratios
    isofox_gene_ids        // channel: [optional]  /path/to/gene_ids
    isofox_tpm_norm        // channel: [optional]  /path/to/tpm_norm

    // Params
    isofox_functions       //  string: [optional]  Isofox functions
    isofox_read_length     //  string: [mandatory] Isofox read length

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources and sort
    // channel: runnable: [ meta, tumor_bam, tumor_bai ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_tumor_rna_bam
        .map { meta, tumor_bam, tumor_bai ->
            return [
                meta,
                Utils.selectCurrentOrExisting(tumor_bam, meta, Constants.INPUT.BAM_RNA_TUMOR),
                Utils.selectCurrentOrExisting(tumor_bai, meta, Constants.INPUT.BAI_RNA_TUMOR),
            ]
        }
        .branch { meta, tumor_bam, tumor_bai ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.ISOFOX_DIR)
            runnable: tumor_bam && !has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_isofox, tumor_bam, tumor_bai ]
    ch_isofox_inputs = ch_inputs_sorted.runnable
        .map { meta, tumor_bam, tumor_bai ->

            def meta_isofox = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorRnaSampleName(meta),
            ]

            return [meta_isofox, tumor_bam, tumor_bai]
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
        known_fusion_data,
        isofox_counts,
        isofox_gc_ratios,
        isofox_gene_ids,
        isofox_tpm_norm,
    )

    ch_versions = ch_versions.mix(ISOFOX.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, isofox_dir ]
    ch_outputs = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ISOFOX.out.isofox_dir, ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    isofox_dir = ch_outputs  // channel: [ meta, isofox_dir ]

    versions   = ch_versions // channel: [ versions.yml ]
}
