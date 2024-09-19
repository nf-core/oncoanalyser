//
// AMBER determines b-allele frequencies at predetermined positions
//

import Constants
import Utils

include { AMBER } from '../../../modules/local/amber/main'

workflow AMBER_PROFILING {
    take:
    // Sample data
    ch_inputs          // channel: [mandatory] [ meta ]
    ch_tumor_bam       // channel: [mandatory] [ meta, bam, bai ]
    ch_normal_bam      // channel: [mandatory] [ meta, bam, bai ]
    ch_donor_bam       // channel: [mandatory] [ meta, bam, bai ]

    // Reference data
    genome_version     // channel: [mandatory] genome version
    heterozygous_sites // channel: [optional]  /path/to/heterozygous_sites
    target_region_bed  // channel: [optional]  /path/to/target_region_bed

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources and sort
    // channel: runnable: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_tumor_bam,
        ch_normal_bam,
        ch_donor_bam,
    )
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai ->
            return [
                meta,
                Utils.selectCurrentOrExisting(tumor_bam, meta, Constants.INPUT.BAM_REDUX_DNA_TUMOR),
                tumor_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_TUMOR),

                Utils.selectCurrentOrExisting(normal_bam, meta, Constants.INPUT.BAM_REDUX_DNA_NORMAL),
                normal_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_NORMAL),

                Utils.selectCurrentOrExisting(donor_bam, meta, Constants.INPUT.BAM_REDUX_DNA_DONOR),
                donor_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_DONOR),
            ]
        }
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.AMBER_DIR)
            runnable: tumor_bam && !has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_amber, tumor_bam, normal_bam, donor_bam, tumor_bai, normal_bai, donor_bai ]
    ch_amber_inputs = ch_inputs_sorted.runnable
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai ->

            def meta_amber = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Utils.getTumorDnaSampleName(meta),
            ]

            if (normal_bam) {
                meta_amber.normal_id = Utils.getNormalDnaSampleName(meta)
            }

            if (donor_bam) {
                meta_amber.donor_id = Utils.getDonorDnaSampleName(meta)
            }

            [meta_amber, tumor_bam, normal_bam, donor_bam, tumor_bai, normal_bai, donor_bai]
        }

    // Run process
    AMBER(
        ch_amber_inputs,
        genome_version,
        heterozygous_sites,
        target_region_bed,
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
