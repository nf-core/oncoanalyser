//
// GRIPSS performs SV filtering.
//

include { GRIPSS_GERMLINE as GERMLINE } from '../../modules/local/gripss/germline/main'
include { GRIPSS_SOMATIC as SOMATIC   } from '../../modules/local/gripss/somatic/main'

workflow GRIPSS_FILTERING {
    take:
        // Sample inputs
        ch_inputs                // channel: [val(meta)]
        ch_gridss                // channel: [val(meta), gridss_vcf]

        // Reference data
        ref_data_genome_fasta    //    file: /path/to/genome_fasta
        ref_data_genome_fai      //    file: /path/to/genome_fai
        ref_data_genome_version  //     val: genome version
        breakend_pon             //    file: /path/to/breakend_pon
        breakpoint_pon           //    file: /path/to/breakpoint_pon
        known_fusions            //    file: /path/to/known_fusions
        repeatmasker_annotations //    file: /path/to/repeatmasker_annotations
        target_regions_bed

        // Params
        run_config

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Select input source
        // channel: [val(meta), gridss_vcf]
        ch_gripss_inputs_source = run_config.stages.gridss ? ch_gridss : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.GRIDSS_VCF)

        // Create inputs and create process-specific meta
        // channel: [val(meta_gripss), gridss_vcf]
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
        ch_germline_out = Channel.empty()
        ch_germline_unfiltered_out = Channel.empty()
        if (run_config.type == Constants.RunType.TUMOR_NORMAL) {

            GERMLINE(
                ch_gripss_inputs,
                ref_data_genome_fasta,
                ref_data_genome_fai,
                ref_data_genome_version,
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
            ref_data_genome_fasta,
            ref_data_genome_fai,
            ref_data_genome_version,
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
        somatic             = ch_somatic_out             // channel: [val(meta), vcf, tbi]
        germline            = ch_germline_out            // channel: [val(meta), vcf, tbi]
        somatic_unfiltered  = ch_somatic_unfiltered_out  // channel: [val(meta), vcf, tbi]
        germline_unfiltered = ch_germline_unfiltered_out // channel: [val(meta), vcf, tbi]

        versions = ch_versions                           // channel: [versions.yml]
}
