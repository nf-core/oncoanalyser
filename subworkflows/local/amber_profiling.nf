//
// AMBER determines b-allele frequencies at predetermined positions
//

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
        ch_inputs_sorted = ch_inputs.branch { meta ->

            //def key = [Constants.SampleType.TUMOR, Constants.SequenceType.DNA]
            //def has_bam = meta.containsKey(key) && meta[key].containsKey(Constants.FileType.BAM)

            ////def has_bam = Utils.hasTumorDnaBam(meta)
            //def existing_input = WorkflowOncoanalyser.getInput(meta, Constants.INPUT.AMBER_DIR)

            existing: Utils.hasExistingInput(meta, Constants.INPUT.AMBER_DIR)
            runnable: Utils.hasTumorDnaBam(meta)
            // NOTE(SW): skip would be included where this entry requires optional AMBER input later on
            //skip: true
        }

        // Create process input channel
        // channel: [ meta_amber, tumor_bam, normal_bam, tumor_bai, normal_bai ]
        ch_amber_inputs = ch_inputs_sorted.runnable
            .map { meta ->

                def tumor_id = Utils.getTumorDnaSampleName(meta)
                def meta_amber = [
                    key: meta.group_id,
                    id: "${meta.group_id}__${tumor_id}",
                    tumor_id: tumor_id,
                ]

                def tumor_bam = Utils.getTumorDnaBam(meta)

                def normal_bam = []
                def normal_bai = []

                def normal_key = [Constants.SampleType.NORMAL, Constants.SequenceType.DNA]
                if (meta.containsKey(normal_key) && meta[normal_key].containsKey(Constants.FileType.BAM)) {

                    meta_amber.normal_id = Utils.getNormalDnaSampleName(meta)
                    normal_bam = Utils.getNormalDnaBam(meta)
                    normal_bai = "${normal_bam}.bai"

                }

                [meta_amber, tumor_bam, normal_bam, "${tumor_bam}.bai", normal_bai]
            }

        // Run process
        AMBER(
            ch_amber_inputs,
            genome_version,
            heterozygous_sites,
        )

        // Set outputs, restoring original meta
        ch_outputs = Channel.empty()
            .mix(
                ch_inputs_sorted.existing,
                WorkflowOncoanalyser.restoreMeta(AMBER.out.amber_dir, ch_inputs),
            )
        ch_versions = ch_versions.mix(AMBER.out.versions)

    emit:
        amber_dir = ch_outputs  // channel: [ meta, amber_dir ]

        versions  = ch_versions // channel: [ versions.yml ]
}
