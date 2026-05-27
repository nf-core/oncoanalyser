//
// AMBER determines b-allele frequencies at predetermined positions
//

include { AMBER } from '../../../modules/local/amber/main'

workflow AMBER_PROFILING {
    take:
    // Sample data
    ch_inputs          // channel: [mandatory] [ meta ]
    ch_tumor_bam       // channel: [mandatory] [ meta, bam, bai ]
    ch_normal_bam      // channel: [mandatory] [ meta, bam, bai ]
    ch_donor_bam       // channel: [mandatory] [ meta, bam, bai ]

    // Reference data
    genome_version       // channel: [mandatory] genome version
    heterozygous_sites   // channel: [optional]  /path/to/heterozygous_sites
    target_regions_bed   // channel: [optional]  /path/to/target_regions_bed
    tumor_min_depth      // integer: [optional]  -tumor_min_depth argument value

    // Params
    sequencing_type      // string:  [mandatory] sequencing type
    purity_estimate_mode // boolean: [mandatory] Set purity estimate mode

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources and sort
    // channel: runnable: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
    // channel: skip: [ meta ]
    ch_inputs_sorted = channels.WorkflowChannels.groupByMeta(
        ch_tumor_bam,
        ch_normal_bam,
        ch_donor_bam,
    )
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai ->
            return [
                meta,
                sample.Inputs.preferUserProvidedInput(tumor_bam, meta, sample.FileKey.BAM_REDUX_DNA_TUMOR),
                sample.Inputs.preferPipelineOutput(tumor_bai, meta, sample.FileKey.BAI_DNA_TUMOR),

                sample.Inputs.preferUserProvidedInput(normal_bam, meta, sample.FileKey.BAM_REDUX_DNA_NORMAL),
                sample.Inputs.preferPipelineOutput(normal_bai, meta, sample.FileKey.BAI_DNA_NORMAL),

                sample.Inputs.preferUserProvidedInput(donor_bam, meta, sample.FileKey.BAM_REDUX_DNA_DONOR),
                sample.Inputs.preferPipelineOutput(donor_bai, meta, sample.FileKey.BAI_DNA_DONOR),
            ]
        }
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai ->
            def has_existing = sample.Inputs.hasExisting(meta, sample.FileKey.AMBER_DIR)

            def runnable_standard = !purity_estimate_mode && tumor_bam && !has_existing

            // TODO(SW): must improve handling through separation of sample information in meta; currently unable to provide ccfDNA AMBER directory in samplesheet
            def runnable_purity_estimate = purity_estimate_mode && normal_bam

            runnable: runnable_standard || runnable_purity_estimate
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_amber, tumor_bam, normal_bam, donor_bam, tumor_bai, normal_bai, donor_bai ]
    ch_amber_inputs = ch_inputs_sorted.runnable
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai ->

            def tumor_id = purity_estimate_mode
                ? sample.Inputs.getTumorDnaSampleNameLongitudinal(meta)
                : sample.Inputs.getTumorDnaSampleNamePrimary(meta)

            def meta_amber = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: tumor_id,
            ]

            if (normal_bam) {
                meta_amber.normal_id = sample.Inputs.getNormalDnaSampleName(meta)
            }

            if (donor_bam) {
                meta_amber.donor_id = sample.Inputs.getDonorDnaSampleName(meta)
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
        sequencing_type,
    )

    ch_versions = ch_versions.mix(AMBER.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, amber_dir ]
    ch_outputs = Channel.empty()
        .mix(
            channels.WorkflowChannels.restoreMeta(AMBER.out.amber_dir, ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    amber_dir = ch_outputs  // channel: [ meta, amber_dir ]

    versions  = ch_versions // channel: [ versions.yml ]
}
