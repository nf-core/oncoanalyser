//
// AMBER determines b-allele frequencies at predetermined positions
//

import Constants
import Inputs

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
    target_regions_bed // channel: [optional]  /path/to/target_regions_bed
    tumor_min_depth    // integer: [optional]  -tumor_min_depth argument value

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
                Inputs.overrideWithExistingInput(tumor_bam, meta, Constants.INPUT.BAM_REDUX_DNA_TUMOR),
                Inputs.fallbackToExistingInput(tumor_bai, meta, Constants.INPUT.BAI_DNA_TUMOR),

                Inputs.overrideWithExistingInput(normal_bam, meta, Constants.INPUT.BAM_REDUX_DNA_NORMAL),
                Inputs.fallbackToExistingInput(normal_bai, meta, Constants.INPUT.BAI_DNA_NORMAL),

                Inputs.overrideWithExistingInput(donor_bam, meta, Constants.INPUT.BAM_REDUX_DNA_DONOR),
                Inputs.fallbackToExistingInput(donor_bai, meta, Constants.INPUT.BAI_DNA_DONOR),
            ]
        }
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai ->
            def has_existing = Inputs.hasExistingInput(meta, Constants.INPUT.AMBER_DIR)


            // TODO(SW): must improve handling through separation of sample information in meta; currently unable to provide ccfDNA AMBER directory in samplesheet
            def longitudinal_sample = Inputs.getTumorDnaSample(meta).containsKey('longitudinal_sample_id')

            runnable: tumor_bam && (!has_existing || longitudinal_sample)


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
                tumor_id: Inputs.getTumorDnaSampleName(meta),
            ]

            if (normal_bam) {
                meta_amber.normal_id = Inputs.getNormalDnaSampleName(meta)
            }

            if (donor_bam) {
                meta_amber.donor_id = Inputs.getDonorDnaSampleName(meta)
            }

            [meta_amber, tumor_bam, normal_bam, donor_bam, tumor_bai, normal_bai, donor_bai]
        }

    // Run process
    AMBER(
        ch_amber_inputs,
        genome_version,
        heterozygous_sites,
        target_regions_bed,
        tumor_min_depth,
    )

    ch_versions = ch_versions.mix(AMBER.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, amber_dir ]
    ch_outputs = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(AMBER.out.amber_dir, ch_inputs),
            PlaceholderChannels.toolDir(ch_inputs_sorted.skip),
        )

    emit:
    amber_dir = ch_outputs  // channel: [ meta, amber_dir ]

    versions  = ch_versions // channel: [ versions.yml ]
}
