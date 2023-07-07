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

        // Params
        run_config               // channel: [mandatory] run configuration

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Select input source
        ch_gripss_inputs_source = run_config.stages.gridss ? ch_gridss : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.GRIDSS_VCF)

        // Create inputs and create process-specific meta
        // channel: [ meta_gripss, gridss_vcf ]
        ch_gripss_inputs = ch_gripss_inputs_source
            .filter { it[0] != Constants.PLACEHOLDER_META }
            .map { meta, gridss_vcf ->

                def meta_gripss = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: Utils.getTumorSampleName(meta, run_config.mode),
                ]

                if (run_config.type == Constants.RunType.TUMOR_NORMAL) {
                    meta_gripss.normal_id = Utils.getNormalWgsSampleName(meta)
                }

                return [meta_gripss, gridss_vcf]
            }

        //
        // MODULE: GRIPSS germline
        //
        // channel: [ meta, gripss_vcf, gripss_tbi ]
        ch_germline_out = Channel.empty()
        ch_germline_unfiltered_out = Channel.empty()
        if (run_config.type == Constants.RunType.TUMOR_NORMAL) {

            GERMLINE(
                ch_gripss_inputs,
                genome_fasta,
                genome_version,
                genome_fai,
                breakend_pon,
                breakpoint_pon,
                known_fusions,
                repeatmasker_annotations,
            )

            ch_versions = ch_versions.mix(GERMLINE.out.versions)
            ch_germline_out = WorkflowOncoanalyser.restoreMeta(GERMLINE.out.vcf, ch_inputs)
            ch_germline_unfiltered_out = WorkflowOncoanalyser.restoreMeta(GERMLINE.out.vcf_unfiltered, ch_inputs)

        }

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
        ch_somatic_out = WorkflowOncoanalyser.restoreMeta(SOMATIC.out.vcf, ch_inputs)
        ch_somatic_unfiltered_out = WorkflowOncoanalyser.restoreMeta(SOMATIC.out.vcf_unfiltered, ch_inputs)

    emit:
        somatic             = ch_somatic_out             // channel: [ meta, gripss_vcf, gripss_tbi ]
        germline            = ch_germline_out            // channel: [ meta, gripss_vcf, gripss_tbi ]
        somatic_unfiltered  = ch_somatic_unfiltered_out  // channel: [ meta, gripss_vcf, gripss_tbi ]
        germline_unfiltered = ch_germline_unfiltered_out // channel: [ meta, gripss_vcf, gripss_tbi ]

        versions = ch_versions                           // channel: [ versions.yml ]
}
