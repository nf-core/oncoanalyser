//
// COBALT determines tumor/normal read ratios for downstream CNV calling
//
import Utils

include { COBALT } from '../../modules/local/cobalt/main'

include { CHANNEL_GROUP_INPUTS } from './channel_group_inputs'

workflow COBALT_PROFILING {
    take:
        // Sample data
        ch_inputs           // channel: [val(meta)]

        // Reference data
        ref_data_gc_profile //    file: /path/to/gc_profile

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Get input meta groups
        CHANNEL_GROUP_INPUTS(
            ch_inputs,
        )

        // Select input sources
        // channel: [meta_cobalt, tumor_bam_wgs, normal_bam_wgs, tumor_bai_wgs, normal_bai_wgs]
        ch_cobalt_inputs = CHANNEL_GROUP_INPUTS.out.wgs_present
            .map { meta ->
                def meta_cobalt = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: Utils.getTumorWgsSampleName(meta),
                    normal_id: Utils.getNormalWgsSampleName(meta),
                ]
                def tumor_bam = Utils.getTumorWgsBam(meta)
                def normal_bam = Utils.getNormalWgsBam(meta)
                return [meta_cobalt, tumor_bam, normal_bam, "${tumor_bam}.bai", "${normal_bam}.bai"]
            }

        // Run process
        COBALT(
            ch_cobalt_inputs,
            ref_data_gc_profile,
        )

        // Set outputs, restoring original meta
        ch_outputs = WorkflowOncoanalyser.restoreMeta(COBALT.out.cobalt_dir, ch_inputs)
        ch_versions = ch_versions.mix(COBALT.out.versions)

    emit:
        cobalt_dir = ch_outputs  // channel: [val(meta), cobalt_dir]

        versions   = ch_versions // channel: [versions.yml]
}
