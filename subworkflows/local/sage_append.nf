//
// SAGE append adds WTS data to an existing SAGE VCF
//

include { SAGE_APPEND as SOMATIC } from '../../modules/local/sage/append/main'
include { SAGE_APPEND as GERMLINE } from '../../modules/local/sage/append/main'

include { CHANNEL_GROUP_INPUTS } from './channel_group_inputs'

workflow SAGE_APPEND {
    take:
        // Sample data
        ch_inputs             // channel: [val(meta)]
        ch_purple_dir         // channel: [val(meta), purple_dir]

        // Reference data
        ref_data_genome_fasta //    file: /path/to/genome_fasta
        ref_data_genome_fai   //    file: /path/to/genome_fai
        ref_data_genome_dict  //    file: /path/to/genome_dict

        // Params
        run

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Get input meta groups
        CHANNEL_GROUP_INPUTS(
            ch_inputs,
        )

        // Select input sources
        // channel: [meta, purple_dir, tumor_wts_bam, tumor_wts_bai]
        ch_sage_append_inputs_source = WorkflowOncoanalyser.groupByMeta(
            run.purple ? ch_purple_dir : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PURPLE_DIR),
            CHANNEL_GROUP_INPUTS.out.wts_present
                .map { meta ->
                    def bam = Utils.getTumorWtsBam(meta)
                    return [meta, bam, "${bam}.bai"]
                },
        )

        // Create inputs and create process-specific meta
        // channel: [sage_meta, purple_somatic_smlv_vcf, tumor_wts_bam]
        ch_sage_somatic_append_inputs = ch_sage_append_inputs_source
            .map { meta, purple_dir, bam, bai ->
                def tumor_id = Utils.getTumorWgsSampleName(meta)
                def purple_smlv_vcf = file(purple_dir).resolve("${tumor_id}.purple.somatic.vcf.gz")

                // Require both somatic smlv from the PURPLE directory
                if (!purple_smlv_vcf.exists()) {
                    return Constants.META_PLACEHOLDER
                }

                def sage_meta = [
                    key: meta.id,
                    id: tumor_id,
                    tumor_wts_id: Utils.getTumorWtsSampleName(meta),
                    wgs_id: tumor_id,
                ]
                return [sage_meta, purple_smlv_vcf, bam, bai]
            }
            .filter { it != Constants.META_PLACEHOLDER }

        // NOTE(SW): revise to reduce repetition

        // channel: [sage_meta, purple_germline_smlv_vcf, tumor_wts_bam]
        ch_sage_germline_append_inputs = ch_sage_append_inputs_source
            .map { meta, purple_dir, bam, bai ->
                def tumor_id = Utils.getTumorWgsSampleName(meta)
                def normal_id = Utils.getNormalWgsSampleName(meta)

                def purple_smlv_vcf = file(purple_dir).resolve("${tumor_id}.purple.germline.vcf.gz")

                // Require both germline smlv from the PURPLE directory
                if (!purple_smlv_vcf.exists()) {
                    return Constants.META_PLACEHOLDER
                }

                def sage_meta = [
                    key: meta.id,
                    id: normal_id,
                    tumor_wts_id: Utils.getTumorWtsSampleName(meta),
                    wgs_id: normal_id,
                ]
                return [sage_meta, purple_smlv_vcf, bam, bai]
            }
            .filter { it != Constants.META_PLACEHOLDER }

        // Run process
        SOMATIC(
            ch_sage_somatic_append_inputs,
            ref_data_genome_fasta,
            ref_data_genome_fai,
            ref_data_genome_dict,
        )

        GERMLINE(
            ch_sage_germline_append_inputs,
            ref_data_genome_fasta,
            ref_data_genome_fai,
            ref_data_genome_dict,
        )

        // Set outputs, restoring original meta
        ch_somatic_vcf = WorkflowOncoanalyser.restoreMeta(SOMATIC.out.vcf, ch_inputs)
        ch_germline_vcf = WorkflowOncoanalyser.restoreMeta(GERMLINE.out.vcf, ch_inputs)
        ch_versions = ch_versions.mix(
          GERMLINE.out.versions,
          SOMATIC.out.versions,
        )

    emit:
        somatic_vcf  = ch_somatic_vcf  // channel: [val(meta), somatic_vcf]
        germline_vcf = ch_germline_vcf // channel: [val(meta), germline_vcf]

        versions     = ch_versions     // channel: [versions.yml]
}
