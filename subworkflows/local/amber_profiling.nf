//
// XXX
//
import Utils

include { AMBER } from '../../modules/local/amber/main'

workflow AMBER_PROFILING {
    take:
        // Sample data
        ch_inputs

        // Reference data
        ref_data_genome_version
        ref_data_heterozygous_sites

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [val(meta_amber), tumor_bam_wgs, normal_bam_wgs, tumor_bai_wgs, normal_bai_wgs]
        ch_amber_inputs = ch_inputs
            .map { meta ->
                def meta_amber = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: Utils.getTumorWgsSampleName(meta),
                    normal_id: Utils.getNormalWgsSampleName(meta),
                ]
                def tumor_bam = Utils.getTumorWgsBam(meta)
                def normal_bam = Utils.getNormalWgsBam(meta)
                [meta_amber, tumor_bam, normal_bam, "${tumor_bam}.bai", "${normal_bam}.bai"]
            }

        // Run process
        AMBER(
            ch_amber_inputs,
            ref_data_genome_version,
            ref_data_heterozygous_sites,
        )

        // Set outputs, restoring original meta
        ch_outputs = WorkflowOncoanalyser.restoreMeta(AMBER.out.amber_dir, ch_inputs)
        ch_versions = ch_versions.mix(AMBER.out.versions)

    emit:
        amber_dir = ch_outputs  // channel: [val(meta), amber_dir]

        versions  = ch_versions // channel: [versions.yml]
}
