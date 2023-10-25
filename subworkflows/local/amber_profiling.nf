//
// AMBER determines b-allele frequencies at predetermined positions
//

import Constants
import Utils

include { AMBER } from '../../modules/local/amber/main'

workflow AMBER_PROFILING {
    take:
        // Sample data
        ch_inputs          // channel: [mandatory] [ meta ]

        // Reference data
        genome_version     // channel: [mandatory] genome version
        heterozygous_sites // channel: [optional]  /path/to/heterozygous_sites

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Sort inputs
        // channel: [ meta ]
        ch_inputs_sorted = ch_inputs
            .branch { meta ->
                def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.AMBER_DIR)
                runnable: Utils.hasTumorDnaBam(meta) && !has_existing
                skip: true
            }

        // Create process input channel
        // channel: [ meta_amber, tumor_bam, normal_bam, tumor_bai, normal_bai ]
        ch_amber_inputs = ch_inputs_sorted.runnable
            .map { meta ->

                def meta_amber = [
                    key: meta.group_id,
                    id: meta.group_id,
                    tumor_id: Utils.getTumorDnaSampleName(meta),
                ]

                def tumor_bam = Utils.getTumorDnaBam(meta)
                def tumor_bai = Utils.getTumorDnaBai(meta)

                def normal_bam = []
                def normal_bai = []

                if (Utils.hasNormalDnaBam(meta)) {

                    meta_amber.normal_id = Utils.getNormalDnaSampleName(meta)
                    normal_bam = Utils.getNormalDnaBam(meta)
                    normal_bai = Utils.getNormalDnaBai(meta)

                }

                [meta_amber, tumor_bam, normal_bam, tumor_bai, normal_bai]
            }

        // Run process
        AMBER(
            ch_amber_inputs,
            genome_version,
            heterozygous_sites,
        )

        ch_versions = ch_versions.mix(AMBER.out.versions)

        // Set outputs, restoring original meta
        // channel: [ meta, amber_dir ]
        ch_outputs = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(AMBER.out.amber_dir, ch_inputs),
                ch_inputs_sorted.skip.map { meta -> [meta, []] },
            )

    emit:
        amber_dir = ch_outputs  // channel: [ meta, amber_dir ]

        versions  = ch_versions // channel: [ versions.yml ]
}
