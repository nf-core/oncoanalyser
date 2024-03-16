//
// GRIPSS performs SV filtering.
//

import Constants
import Utils

include { GRIPSS_GERMLINE as GERMLINE } from '../../../modules/local/gripss/germline/main'
include { GRIPSS_SOMATIC as SOMATIC   } from '../../../modules/local/gripss/somatic/main'

workflow GRIPSS_FILTERING {
    take:
        // Sample inputs
        ch_inputs                // channel: [mandatory] [ meta ]
        ch_gridss                // channel: [mandatory] [ meta, gridss_vcf ]

        // Reference data
        genome_fasta             // channel: [mandatory] /path/to/genome_fasta
        genome_version           // channel: [mandatory] genome version
        genome_fai               // channel: [mandatory] /path/to/genome_fai
        breakend_pon             // channel: [mandatory] /path/to/breakend_pon
        breakpoint_pon           // channel: [mandatory] /path/to/breakpoint_pon
        known_fusions            // channel: [mandatory] /path/to/known_fusions
        repeatmasker_annotations // channel: [mandatory] /path/to/repeatmasker_annotations
        target_region_bed        // channel: [optional]  /path/to/target_region_bed

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Select input sources and sort
        // channel: runnable: [ meta, gridss_vcf ]
        // channel: skip: [ meta ]
        ch_inputs_sorted = ch_gridss
            .map { meta, gridss_vcf ->
                return [
                    meta,
                    Utils.selectCurrentOrExisting(gridss_vcf, meta, Constants.INPUT.GRIDSS_VCF),
                ]
            }
            .branch { meta, gridss_vcf ->
                runnable: gridss_vcf
                skip: true
                    return meta
            }

        //
        // MODULE: GRIPSS germline
        //
        // Select inputs that are eligible to run
        // channel: runnable: [ meta, gridss_vcf ]
        // channel: skip: [ meta ]
        ch_inputs_germline_sorted = ch_inputs_sorted.runnable
            .branch { meta, gridss_vcf ->
                def has_tumor_normal = Utils.hasTumorDna(meta) && Utils.hasNormalDna(meta)
                def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.GRIPSS_VCF_NORMAL)

                runnable: has_tumor_normal && !has_existing
                skip: true
                    return meta
            }

        // Create process input channel
        // channel: [ meta_gripss, gridss_vcf ]
        ch_gripss_germline_inputs = ch_inputs_germline_sorted.runnable
            .map { meta, gridss_vcf ->

                def meta_gripss = [
                    key: meta.group_id,
                    id: meta.group_id,
                    tumor_id: Utils.getTumorDnaSampleName(meta),
                    normal_id: Utils.getNormalDnaSampleName(meta),
                ]

                return [meta_gripss, gridss_vcf]
            }

        // Run process
        GERMLINE(
            ch_gripss_germline_inputs,
            genome_fasta,
            genome_version,
            genome_fai,
            breakend_pon,
            breakpoint_pon,
            known_fusions,
            repeatmasker_annotations,
        )

        ch_versions = ch_versions.mix(GERMLINE.out.versions)

        //
        // MODULE: GRIPSS somatic
        //
        // Select inputs that are eligible to run
        // channel: runnable: [ meta, gridss_vcf ]
        // channel: skip: [ meta ]
        ch_inputs_somatic_sorted = ch_inputs_sorted.runnable
            .branch { meta, gridss_vcf ->
                def has_tumor = Utils.hasTumorDna(meta)
                def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.GRIPSS_VCF_TUMOR)

                runnable: has_tumor && !has_existing
                skip: true
                    return meta
            }

        // Create process input channel
        // channel: [ meta_gripss, gridss_vcf ]
        ch_gripss_somatic_inputs = ch_inputs_somatic_sorted.runnable
            .map { meta, gridss_vcf ->

                def meta_gripss = [
                    key: meta.group_id,
                    id: meta.group_id,
                    tumor_id: Utils.getTumorDnaSampleName(meta),
                ]

                if (Utils.hasNormalDna(meta)) {
                    meta_gripss.normal_id = Utils.getNormalDnaSampleName(meta)
                }

                return [meta_gripss, gridss_vcf]
            }

        // Run process
        SOMATIC(
            ch_gripss_somatic_inputs,
            genome_fasta,
            genome_version,
            genome_fai,
            breakend_pon,
            breakpoint_pon,
            known_fusions,
            repeatmasker_annotations,
            target_region_bed,
        )

        ch_versions = ch_versions.mix(SOMATIC.out.versions)

        // Set outputs, restoring original meta
        // channel: [ meta, gripss_vcf, gripss_tbi ]
        ch_somatic_out = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(SOMATIC.out.vcf, ch_inputs),
                ch_inputs_somatic_sorted.skip.map { meta -> [meta, [], []] },
                ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
            )

        ch_somatic_unfiltered_out = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(SOMATIC.out.vcf_unfiltered, ch_inputs),
                ch_inputs_somatic_sorted.skip.map { meta -> [meta, [], []] },
                ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
            )

        ch_germline_out = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(GERMLINE.out.vcf, ch_inputs),
                ch_inputs_germline_sorted.skip.map { meta -> [meta, [], []] },
                ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
            )

        ch_germline_unfiltered_out = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(GERMLINE.out.vcf_unfiltered, ch_inputs),
                ch_inputs_germline_sorted.skip.map { meta -> [meta, [], []] },
                ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
            )

    emit:
        somatic             = ch_somatic_out             // channel: [ meta, gripss_vcf, gripss_tbi ]
        germline            = ch_germline_out            // channel: [ meta, gripss_vcf, gripss_tbi ]
        somatic_unfiltered  = ch_somatic_unfiltered_out  // channel: [ meta, gripss_vcf, gripss_tbi ]
        germline_unfiltered = ch_germline_unfiltered_out // channel: [ meta, gripss_vcf, gripss_tbi ]

        versions = ch_versions                           // channel: [ versions.yml ]
}
