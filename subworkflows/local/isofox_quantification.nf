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
        run_config             // channel: [mandatory] run configuration

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Create inputs and create process-specific meta
        // channel: [ meta_isofox, tumor_bam_rna ]
        if (run_config.stages.isofox) {
            ch_isofox_inputs = ch_inputs
                .map { meta ->
                    def bam = Utils.getTumorRnaBam(meta)
                    def meta_isofox = [key: meta.id, id: Utils.getTumorRnaSampleName(meta)]
                    return [meta_isofox, bam, "${bam}.bai"]
                }
        } else {
            ch_isofox_inputs = WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.ISOFOX_DIR, type: 'optional')
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

        // Set outputs, restoring original meta
        ch_outputs = WorkflowOncoanalyser.restoreMeta(ISOFOX.out.isofox_dir, ch_inputs)
        ch_versions = ch_versions.mix(ISOFOX.out.versions)

    emit:
        isofox_dir = ch_outputs // channel: [ meta, isofox_dir ]

        versions  = ch_versions // channel: [ versions.yml ]
}
