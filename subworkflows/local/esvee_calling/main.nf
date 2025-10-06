//
// ESVEE detects structural variants, and reports breakends and breakpoints.
//

import Constants
import Utils

import java.nio.channels.Channel

include { ESVEE } from '../../../modules/local/esvee/main'

workflow ESVEE_CALLING {
    take:
    ch_inputs                // channel: [mandatory] [ meta ]
    ch_tumor_bam             // channel: [mandatory] [ meta, bam, bai ]
    ch_normal_bam            // channel: [mandatory] [ meta, bam, bai ]
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

    main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Select input sources and sort
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_tumor_bam,
        ch_normal_bam,
    )
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            return [
                meta,
                Utils.selectCurrentOrExisting(tumor_bam, meta, Constants.INPUT.BAM_REDUX_DNA_TUMOR),
                tumor_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_TUMOR),
                Utils.selectCurrentOrExisting(normal_bam, meta, Constants.INPUT.BAM_REDUX_DNA_NORMAL),
                normal_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_NORMAL),
            ]
        }
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.ESVEE_VCF_TUMOR)

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
                tumor_id: Utils.getTumorDnaSampleName(meta),
            ]

            if (normal_bam) {
                meta_esvee.normal_id = Utils.getNormalDnaSampleName(meta)
            }

            return [meta_esvee, tumor_bam, tumor_bai, normal_bam, normal_bai]
        }

    // Run ESVEE process
    ESVEE(
        ch_esvee_inputs,
        genome_fasta,
        genome_fai,
        genome_dict,
        genome_img,
        genome_version,
        pon_breakends,
        pon_breakpoints,
        decoy_sequences_image,
        known_fusions,
        repeatmasker_annotations,
        unmap_regions,
        target_region_bed
    )

    ch_versions = ch_versions.mix(ESVEE.out.versions)

    // Set outputs, restoring original meta
    ch_somatic_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ESVEE.out.somatic_vcf, ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, [], []] }
        )

    ch_germline_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ESVEE.out.germline_vcf, ch_inputs),
            ch_inputs_sorted.runnable_to.map { meta, tumor_bam, tumor_bai -> [meta, [], []] },
            ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
        )

    ch_unfiltered_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ESVEE.out.unfiltered_vcf, ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, [], []] }
        )

    emit:
    somatic_vcf    = ch_somatic_out    // channel: [ meta, vcf, tbi ]
    germline_vcf   = ch_germline_out   // channel: [ meta, vcf, tbi ]
    unfiltered_vcf = ch_unfiltered_out // channel: [ meta, vcf, tbi ]

    versions       = ch_versions       // channel: [ versions.yml ]
}
