//
// XXX
//
import Utils

include { ISOFOX } from '../../modules/local/isofox/main'

workflow ISOFOX_QUANTIFICATION {
    take:
        // Sample data
        ch_inputs

        // Reference data
        ref_data_genome_fasta
        ref_data_genome_fai
        ref_data_genome_version
        ref_data_ensembl_data_resources
        ref_data_isofox_counts
        ref_data_isofox_gc_ratios

        // Parameters
        isofox_functions
        //use_isofox_exp_counts_cache

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Create inputs and create process-specific meta
        // channel: [meta_isofox, tumor_bam_wts]
        ch_isofox_inputs = ch_inputs
            .map { meta ->
                def bam = Utils.getTumorWtsBam(meta)
                def meta_isofox = [key: meta.id, id: Utils.getTumorWtsSampleName(meta)]
                return [meta_isofox, bam, "${bam}.bai"]
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

        //isofox_counts = params.use_isofox_exp_counts_cache ? ref_data_isofox_counts : []
        isofox_counts = ref_data_isofox_counts

        // Run process
        ISOFOX(
            ch_isofox_inputs,
            isofox_functions,
            ref_data_genome_fasta,
            ref_data_genome_fai,
            ref_data_genome_version,
            ref_data_ensembl_data_resources,
            isofox_counts,
            ref_data_isofox_gc_ratios,
        )

        // Set outputs, restoring original meta
        ch_outputs = WorkflowOncoanalyser.restoreMeta(ISOFOX.out.isofox_dir, ch_inputs)
        ch_versions = ch_versions.mix(ISOFOX.out.versions)

    emit:
        isofox_dir = ch_outputs // channel: [val(meta), isofox_dir]

        versions  = ch_versions // channel: [versions.yml]
}
