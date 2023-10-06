//
// GRIPSS performs SV filtering.
//

include { GRIPSS_GERMLINE as GERMLINE } from '../../modules/local/gripss/germline/main'
include { GRIPSS_SOMATIC as SOMATIC   } from '../../modules/local/gripss/somatic/main'

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
        target_regions_bed       // channel: [optional]  /path/to/target_regions_bed

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Select input sources and sort
        ch_inputs_sorted = ch_gridss
            .map { meta, gridss_vcf ->
                if (Utils.hasExistingInput(meta, Constants.INPUT.GRIDSS_VCF)) {
                    return [meta, Utils.getInput(meta, Constants.INPUT.GRIDSS_VCF)]
                } else {
                    return [meta, gridss_vcf]
                }
            }
            .branch {
                // Runnable status conditioned on presence of gridss_vcf (it[1])
                runnable: it[1]
                skip: true
            }

        // Create general process input channel
        // channel: [ meta_gripss, gridss_vcf ]
        ch_gripss_inputs = ch_inputs_sorted.runnable
            .map { meta, gridss_vcf ->

                // NOTE(SW): germline only is not currently supported
                assert Utils.hasTumorDnaBam(meta)

                def tumor_id = Utils.getTumorDnaSampleName(meta)
                def meta_gripss = [
                    key: meta.group_id,
                    id: "${meta.group_id}__${tumor_id}",
                    tumor_id: tumor_id,
                ]

                if (Utils.hasNormalDnaBam(meta)) {
                    meta_gripss.normal_id = Utils.getNormalDnaSampleName(meta)
                    meta_gripss.sample_type = 'tumor_normal'
                } else if (Utils.hasTumorDnaBam(meta)) {
                    meta_gripss.sample_type = 'tumor_only'
                } else {
                    assert false
                }

                return [meta_gripss, gridss_vcf]
            }

        //
        // MODULE: GRIPSS germline
        //
        // Create germline input channel
        // channel: [ meta, gripss_vcf ]
        ch_gripss_germline_inputs = ch_gripss_inputs
            .filter {
                def meta_gridss = it[0]
                return meta_gridss.sample_type == 'tumor_normal'
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
        SOMATIC(
            ch_gripss_inputs,
            genome_fasta,
            genome_version,
            genome_fai,
            breakend_pon,
            breakpoint_pon,
            known_fusions,
            repeatmasker_annotations,
            target_regions_bed,
        )

        ch_versions = ch_versions.mix(SOMATIC.out.versions)

        // Set outputs, restoring original meta
        // NOTE(SW): look to reduce repetitiveness here
        // channel: [ meta, gripss_vcf, gripss_tbi ]
        ch_somatic_out = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(SOMATIC.out.vcf, ch_inputs),
                ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
            )

        ch_somatic_unfiltered_out = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(SOMATIC.out.vcf_unfiltered, ch_inputs),
                ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
            )

        ch_germline_out = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(GERMLINE.out.vcf, ch_inputs),
                ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
            )

        ch_germline_unfiltered_out = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(GERMLINE.out.vcf_unfiltered, ch_inputs),
                ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
            )

    emit:
        somatic             = ch_somatic_out             // channel: [ meta, gripss_vcf, gripss_tbi ]
        germline            = ch_germline_out            // channel: [ meta, gripss_vcf, gripss_tbi ]
        somatic_unfiltered  = ch_somatic_unfiltered_out  // channel: [ meta, gripss_vcf, gripss_tbi ]
        germline_unfiltered = ch_germline_unfiltered_out // channel: [ meta, gripss_vcf, gripss_tbi ]

        versions = ch_versions                           // channel: [ versions.yml ]
}
