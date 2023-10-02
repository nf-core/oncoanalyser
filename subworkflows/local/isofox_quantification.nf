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

        // Params
        isofox_functions       //  string: [optional] isofox functions
        //use_isofox_exp_counts_cache
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

        // Set Isofox cache files
        // NOTE(SW): the Isofox expected count file is read length dependent so required users to explicitly use expect
        // counts generated for 151 bp reads that is available in the HMF reference bundle. When not specifying an
        // expected count file, Isofox will automatically create one for the computed read length. However, doing so
        // greatly increases runtime.
        // NOTE(SW): consider alternative approaches for using the expected count file e.g. generate once at runtime,
        // then use for all samples; generate all possible read lengths outside of pipeline and store on a remote for
        // retrieval at runtime (requires inference of read length)

        // TODO(SW): this must be improved to allow users to set input file, use cache, or generate at runtime;
        // currently does not update functions
        // NOTE(SW): forcing use of cache for now since this feature is incomplete

        //isofox_counts = params.use_isofox_exp_counts_cache ? isofox_counts : []
        isofox_counts = isofox_counts

        // Run process
        ISOFOX(
            ch_isofox_inputs,
            isofox_functions,
            genome_fasta,
            genome_version,
            genome_fai,
            ensembl_data_resources,
            isofox_counts,
            isofox_gc_ratios,
        )

        // Set outputs, restoring original meta
        ch_outputs = WorkflowOncoanalyser.restoreMeta(ISOFOX.out.isofox_dir, ch_inputs)
        ch_versions = ch_versions.mix(ISOFOX.out.versions)

    emit:
        isofox_dir = ch_outputs // channel: [ meta, isofox_dir ]

        versions  = ch_versions // channel: [ versions.yml ]
}
