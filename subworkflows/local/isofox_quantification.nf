//
// Isofox estimates transcript abundance, detects novel SJs, and identifies fusion events
//

import Utils

include { ISOFOX } from '../../modules/local/isofox/main'

workflow ISOFOX_QUANTIFICATION {
    take:
        // Sample data
        ch_inputs              // channel: [mandatory] [ meta ]

        // Reference data
        genome_fasta           // channel: [mandatory] /path/to/genome_fasta
        genome_version         // channel: [mandatory] genome version
        genome_fai             // channel: [mandatory] /path/to/genome_fai
        ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/
        isofox_counts          // channel: [mandatory] /path/to/isofox_counts
        isofox_gc_ratios       // channel: [mandatory] /path/to/isofox_gc_ratios
        isofox_gene_ids        // channel: [optional]  /path/to/gene_ids
        isofox_tpm_norm        // channel: [optional]  /path/to/tpm_norm

        // Params
        isofox_functions       //  string: [optional]  isofox functions
        isofox_read_length     //  string: [mandatory] isofox_read_length

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Sort inputs
        // channel: [ meta ]
        ch_inputs_sorted = ch_inputs.branch { meta ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.ISOFOX_DIR)
            runnable: Utils.hasTumorRnaBam(meta) && !has_existing
            skip: true
        }

        // Create process input channel
        // channel: [ meta_isofox, tumor_bam_rna ]
        ch_isofox_inputs = ch_inputs_sorted.runnable
            .map { meta ->

                def tumor_id = Utils.getTumorRnaSampleName(meta)
                def meta_isofox = [
                    key: meta.group_id,
                    id: "${meta.group_id}__${tumor_id}",
                ]

                return [meta_isofox, Utils.getTumorRnaBam(meta), Utils.getTumorRnaBai(meta)]
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
        isofox_dir = ch_outputs // channel: [ meta, isofox_dir ]

        versions  = ch_versions // channel: [ versions.yml ]
}
