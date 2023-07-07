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

        // Params
        run_config         // channel: [mandatory] run configuration

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [ meta_amber, tumor_bam, normal_bam, tumor_bai, normal_bai ]
        ch_amber_inputs = ch_inputs
            .map { meta ->
                def meta_amber = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: Utils.getTumorSampleName(meta, run_config.mode),
                ]

                def tumor_bam = Utils.getTumorBam(meta, run_config.mode)

                def normal_bam = []
                def normal_bai = []

                if (run_config.type == Constants.RunType.TUMOR_NORMAL) {

                    assert [Constants.RunMode.WGS, Constants.RunMode.WGTS].contains(run_config.mode)

                    meta_amber.normal_id = Utils.getNormalWgsSampleName(meta)
                    normal_bam = Utils.getNormalWgsBam(meta)
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
        ch_outputs = WorkflowOncoanalyser.restoreMeta(AMBER.out.amber_dir, ch_inputs)
        ch_versions = ch_versions.mix(AMBER.out.versions)

    emit:
        amber_dir = ch_outputs  // channel: [ meta, amber_dir ]

        versions  = ch_versions // channel: [ versions.yml ]
}
