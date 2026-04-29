//
// COBALT calculates read ratios between tumor and normal samples
//

include { COBALT } from '../../../modules/local/cobalt/run/main'

workflow COBALT_PROFILING {
    take:
    // Sample data
    ch_inputs                   // channel: [mandatory] [ meta ]
    ch_tumor_bam                // channel: [mandatory] [ meta, bam, bai ]
    ch_normal_bam               // channel: [mandatory] [ meta, bam, bai ]

    // Reference data
    genome_version              // channel: [mandatory] genome version
    gc_profile                  // channel: [mandatory] /path/to/gc_profile
    diploid_bed                 // channel: [optional]  /path/to/diploid_bed
    target_region_normalisation // channel: [optional]  /path/to/target_region_normalisation
    targeted_mode               // boolean: [mandatory] Set targeted mode
    purity_estimate_mode        // boolean: [mandatory] Set purity estimate mode

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources and sort
    // NOTE(SW): germline mode is not currently supported
    // channel: runnable: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
    // channel: skip: [ meta ]
    ch_inputs_sorted = channels.WorkflowChannels.groupByMeta(
        ch_tumor_bam,
        ch_normal_bam,
    )
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            return [
                meta,
                sample.Inputs.preferUserProvidedInput(tumor_bam, meta, sample.FileKey.BAM_REDUX_DNA_TUMOR),
                sample.Inputs.preferPipelineOutput(tumor_bai, meta, sample.FileKey.BAI_DNA_TUMOR),
                sample.Inputs.preferUserProvidedInput(normal_bam, meta, sample.FileKey.BAM_REDUX_DNA_NORMAL),
                sample.Inputs.preferPipelineOutput(normal_bai, meta, sample.FileKey.BAI_DNA_NORMAL),
            ]
        }
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            def has_existing = sample.Inputs.hasExisting(meta, sample.FileKey.COBALT_DIR)
            runnable_tn: tumor_bam && normal_bam && !has_existing
            runnable_to: tumor_bam && !has_existing
            skip: true
                return meta
        }

    // First set diploid BED input for tumor/normal and tumor only samples
    // NOTE(SW): since the diploid BED is provided as a channel, I seem to be only able to include via channel ops
    // channel: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai, diploid_bed ]
    ch_inputs_runnable = Channel.empty()
        .mix(
            ch_inputs_sorted.runnable_tn.map { [*it, []] },
            ch_inputs_sorted.runnable_to.combine(diploid_bed),
        )

    // Create process input channel
    // channel: sample_data: [ meta_cobalt, tumor_bam, normal_bam, tumor_bai, normal_bai ]
    // channel: diploid_bed: [ diploid_bed ]
    ch_cobalt_inputs = ch_inputs_runnable
        .multiMap { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, diploid_bed ->

            def tumor_id = purity_estimate_mode
               ? sample.Inputs.getTumorDnaSampleName(meta, 'longitudinal')
               : sample.Inputs.getTumorDnaSampleName(meta, 'primary')

            def meta_cobalt = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: tumor_id,
            ]

            if (normal_bam) {
                meta_cobalt.normal_id = sample.Inputs.getNormalDnaSampleName(meta)
            }

            sample_data: [meta_cobalt, tumor_bam, normal_bam, tumor_bai, normal_bai]
            diploid_bed: diploid_bed
        }

    // Run process
    COBALT(
        ch_cobalt_inputs.sample_data,
        genome_version,
        gc_profile,
        ch_cobalt_inputs.diploid_bed,
        target_region_normalisation,
        targeted_mode,
    )

    ch_versions = ch_versions.mix(COBALT.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, cobalt_dir ]
    ch_outputs = Channel.empty()
        .mix(
            channels.WorkflowChannels.restoreMeta(COBALT.out.cobalt_dir, ch_inputs),
            channels.PlaceholderChannels.toolDir(ch_inputs_sorted.skip),
        )

    emit:
    cobalt_dir = ch_outputs  // channel: [ meta, cobalt_dir ]

    versions   = ch_versions // channel: [ versions.yml ]
}
