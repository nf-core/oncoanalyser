//
// SAGE append adds WTS data to an existing SAGE VCF
//

import Constants

include { SAGE_APPEND as SOMATIC } from '../../modules/local/sage/append/main'
include { SAGE_APPEND as GERMLINE } from '../../modules/local/sage/append/main'

workflow SAGE_APPEND {
    take:
        // Sample data
        ch_inputs      // channel: [mandatory] [ meta ]
        ch_purple_dir  // channel: [mandatory] [ meta, purple_dir ]

        // Reference data
        genome_fasta   // channel: [mandatory] /path/to/genome_fasta
        genome_version // channel: [mandatory] genome version
        genome_fai     // channel: [mandatory] /path/to/genome_fai
        genome_dict    // channel: [mandatory] /path/to/genome_dict

        // Params
        run_config     // channel: [mandatory] run configuration

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [ meta, purple_dir, tumor_wts_bam, tumor_wts_bai ]
        ch_sage_append_inputs_source = WorkflowOncoanalyser.groupByMeta(
            run_config.stages.purple ? ch_purple_dir : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PURPLE_DIR),
            ch_inputs
                .map { meta ->
                    def bam = Utils.getTumorWtsBam(meta)
                    return [meta, bam, "${bam}.bai"]
                },
        )

        //
        // MODULE: SAGE append germline
        //
        // channel: [ meta, sage_append_vcf ]
        ch_germline_vcf = Channel.empty()
        if (run_config.type == Constants.RunType.TUMOR_NORMAL) {

            // channel: [ sage_meta, purple_germline_smlv_vcf, tumor_wts_bam ]
            ch_sage_germline_append_inputs = ch_sage_append_inputs_source
                .map { meta, purple_dir, bam, bai ->
                    def tumor_id = Utils.getTumorWgsSampleName(meta)
                    def normal_id = Utils.getNormalWgsSampleName(meta)

                    def purple_smlv_vcf = file(purple_dir).resolve("${tumor_id}.purple.germline.vcf.gz")

                    // Require both germline smlv from the PURPLE directory
                    if (!purple_smlv_vcf.exists()) {
                        return Constants.PLACEHOLDER_META
                    }

                    def sage_meta = [
                        key: meta.id,
                        id: meta.id,
                        tumor_wts_id: Utils.getTumorWtsSampleName(meta),
                        wgs_id: normal_id,
                    ]
                    return [sage_meta, purple_smlv_vcf, bam, bai]
                }
                .filter { it != Constants.PLACEHOLDER_META }

            GERMLINE(
                ch_sage_germline_append_inputs,
                genome_fasta,
                genome_version,
                genome_fai,
                genome_dict,
            )

            // Set outputs, restoring original meta
            ch_germline_vcf = WorkflowOncoanalyser.restoreMeta(GERMLINE.out.vcf, ch_inputs)
            ch_versions = ch_versions.mix(GERMLINE.out.versions)
        }

        //
        // MODULE: SAGE append germline
        //
        // NOTE(SW): revise to reduce repetition
        // Create inputs and create process-specific meta
        // channel: [ sage_meta, purple_somatic_smlv_vcf, tumor_wts_bam ]
        ch_sage_somatic_append_inputs = ch_sage_append_inputs_source
            .map { meta, purple_dir, bam, bai ->
                def tumor_id = Utils.getTumorWgsSampleName(meta)
                def purple_smlv_vcf = file(purple_dir).resolve("${tumor_id}.purple.somatic.vcf.gz")

                // Require both somatic smlv from the PURPLE directory
                if (!purple_smlv_vcf.exists()) {
                    return Constants.PLACEHOLDER_META
                }

                def sage_meta = [
                    key: meta.id,
                    id: meta.id,
                    tumor_wts_id: Utils.getTumorWtsSampleName(meta),
                    wgs_id: tumor_id,
                ]
                return [sage_meta, purple_smlv_vcf, bam, bai]
            }
            .filter { it != Constants.PLACEHOLDER_META }

        SOMATIC(
            ch_sage_somatic_append_inputs,
            genome_fasta,
            genome_version,
            genome_fai,
            genome_dict,
        )

        // Set outputs, restoring original meta
        ch_somatic_vcf = WorkflowOncoanalyser.restoreMeta(SOMATIC.out.vcf, ch_inputs)
        ch_versions = ch_versions.mix(SOMATIC.out.versions)

    emit:
        somatic_vcf  = ch_somatic_vcf  // channel: [ meta, sage_append_vcf ]
        germline_vcf = ch_germline_vcf // channel: [ meta, sage_append_vcf ]

        versions     = ch_versions     // channel: [ versions.yml ]
}
