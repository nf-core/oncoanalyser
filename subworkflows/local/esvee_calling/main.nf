//
// ESVEE detects structural variants, and reports breakends and breakpoints.
//

include { ESVEE } from '../../../modules/local/esvee/main'

workflow ESVEE_CALLING {
    take:

    // Sample data
    ch_inputs                // channel: [mandatory] [ meta ]
    ch_tumor_bam             // channel: [mandatory] [ meta, bam, bai ]
    ch_normal_bam            // channel: [mandatory] [ meta, bam, bai ]

    // Reference data
    genome_fasta             // channel: [mandatory] /path/to/genome_fasta
    genome_version           // channel: [mandatory] genome version
    genome_fai               // channel: [mandatory] /path/to/genome_fai
    genome_dict              // channel: [mandatory] /path/to/genome_dict
    genome_img               // channel: [optional]  /path/to/genome_img
    known_fusions            // channel: [mandatory] /path/to/known_fusions
    pon_breakends            // channel: [mandatory] /path/to/pon_sgl
    pon_breakpoints          // channel: [mandatory] /path/to/pon_sv
    decoy_sequences_image    // channel: [mandatory] /path/to/decoy_sequences_image
    repeatmasker_annotations // channel: [mandatory] /path/to/repeatmasker_annotations
    unmap_regions            // channel: [mandatory] /path/to/unmap_regions
    target_region_bed        // channel: [optional]  /path/to/target_region_bed

    // Params
    sequencing_type          // string:  [mandatory] sequencing type

    main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Select input sources and sort
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
            def has_existing = sample.Inputs.hasExisting(meta, sample.FileKey.ESVEE_DIR)

            runnable_tn: tumor_bam && normal_bam && !has_existing
            runnable_to: tumor_bam && !has_existing
                return [meta, tumor_bam, tumor_bai]
            skip: true
                return meta
        }

    // Create process input channel
    ch_esvee_inputs = Channel.empty()
        .mix(
            ch_inputs_sorted.runnable_tn,
            ch_inputs_sorted.runnable_to.map { [*it, [], []] },
        )
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->

            def meta_esvee = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: sample.Inputs.getTumorDnaSampleName(meta),
            ]

            if (normal_bam) {
                meta_esvee.normal_id = sample.Inputs.getNormalDnaSampleName(meta)
            }

            return [meta_esvee, tumor_bam, tumor_bai, normal_bam, normal_bai]
        }

    // Run ESVEE process
    ESVEE(
        ch_esvee_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        genome_img,
        pon_breakends,
        pon_breakpoints,
        decoy_sequences_image,
        known_fusions,
        repeatmasker_annotations,
        unmap_regions,
        target_region_bed,
        sequencing_type,
    )

    ch_versions = ch_versions.mix(ESVEE.out.versions)

    // Set outputs, restoring original meta
    ch_esvee_out = Channel.empty()
        .mix(
            channels.WorkflowChannels.restoreMeta(ESVEE.out.esvee_dir, ch_inputs),
            channels.PlaceholderChannels.toolDir(ch_inputs_sorted.skip),
        )

    emit:
    esvee_dir      = ch_esvee_out      // channel: [ meta, dir ]

    versions       = ch_versions       // channel: [ versions.yml ]
}
