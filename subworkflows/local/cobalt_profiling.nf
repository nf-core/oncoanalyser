//
// COBALT calculates read ratios between tumor and normal samples
//

import Constants
import Utils

include { COBALT } from '../../modules/local/cobalt/main'

workflow COBALT_PROFILING {
    take:
        // Sample data
        ch_inputs                   // channel: [mandatory] [ meta ]

        // Reference data
        gc_profile                  // channel: [mandatory] /path/to/gc_profile
        diploid_bed                 // channel: [optional]  /path/to/diploid_bed
        target_region_normalisation // channel: [optional]  /path/to/target_region_normalisation

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Sort inputs
        // NOTE(SW): germline mode is not currently supported
        // channel: [ meta ]
        ch_inputs_sorted = ch_inputs.branch { meta ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.COBALT_DIR)
            runnable_tn: Utils.hasTumorDnaBam(meta) && Utils.hasNormalDnaBam(meta) && !has_existing
            runnable_to: Utils.hasTumorDnaBam(meta) && !has_existing
            skip: true
        }

        // First set diploid BED input for tumor/normal and tumor only samples
        // NOTE(SW): since the diploid BED is provided as a channel, I seem to be only able to include via channel ops
        // channel: [ meta, diploid_bed ]
        ch_inputs_runnable = Channel.empty()
            .mix(
                ch_inputs_sorted.runnable_tn.map { meta -> [meta, []] },
                ch_inputs_sorted.runnable_to.combine(diploid_bed),
            )

        // Select input sources
        // channel: [ meta_cobalt, tumor_bam, normal_bam, tumor_bai, normal_bai, diploid_bed ]
        ch_cobalt_inputs = ch_inputs_runnable
            .map { meta, diploid_bed ->

                def tumor_id = Utils.getTumorDnaSampleName(meta)
                def meta_cobalt = [
                    key: meta.group_id,
                    id: "${meta.group_id}__${tumor_id}",
                    tumor_id: tumor_id,
                ]

                def tumor_bam = Utils.getTumorDnaBam(meta)
                def tumor_bai = Utils.getTumorDnaBai(meta)

                def normal_bam = []
                def normal_bai = []

                if (Utils.hasNormalDnaBam(meta)) {
                    meta_cobalt.normal_id = Utils.getNormalDnaSampleName(meta)
                    normal_bam = Utils.getNormalDnaBam(meta)
                    normal_bai = Utils.getNormalDnaBai(meta)
                }

                return [meta_cobalt, tumor_bam, normal_bam, tumor_bai, normal_bai, diploid_bed]
            }

        // Run process
        COBALT(
            ch_cobalt_inputs,
            gc_profile,
            target_region_normalisation,
        )

        ch_versions = ch_versions.mix(COBALT.out.versions)

        // Set outputs, restoring original meta
        // channel: [ meta, cobalt_dir ]
        ch_outputs = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(COBALT.out.cobalt_dir, ch_inputs),
                ch_inputs_sorted.skip.map { meta -> [meta, []] },
            )

    emit:
        cobalt_dir = ch_outputs  // channel: [ meta, cobalt_dir ]

        versions   = ch_versions // channel: [ versions.yml ]
}
