//
// SAGE is a precise and highly sensitive somatic SNV, MNV and small INDEL caller.
//

include { SAGE_GERMLINE as GERMLINE } from '../../modules/local/sage/germline/main'
include { SAGE_SOMATIC as SOMATIC   } from '../../modules/local/sage/somatic/main'

include { CHANNEL_GROUP_INPUTS } from './channel_group_inputs'

workflow SAGE_CALLING {
    take:
        ch_inputs                             // channel: [val(meta)]
        ref_data_genome_fasta                 //    file: /path/to/genome_fasta
        ref_data_genome_fai                   //    file: /path/to/genome_fai
        ref_data_genome_dict                  //    file: /path/to/genome_dict
        ref_data_genome_version               //     val: genome version
        ref_data_sage_known_hotspots_germline //    file: /path/to/sage_known_hotspots_germline
        ref_data_sage_known_hotspots_somatic  //    file: /path/to/sage_known_hotspots_somatic
        ref_data_sage_actionable_panel        //    file: /path/to/sage_actionable_panel
        ref_data_sage_coverage_panel          //    file: /path/to/sage_coverage_panel
        ref_data_sage_highconf_regions        //    file: /path/to/sage_highconf_regions
        ref_data_sage_pon                     //    file: /path/to/sage_pon
        ref_data_segment_mappability          //    file: /path/to/segment_mappability
        ref_data_driver_gene_panel            //    file: /path/to/driver_gene_panel
        ref_data_ensembl_data_resources       //    file: /path/to/ensembl_data_resources/

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Get input meta groups
        CHANNEL_GROUP_INPUTS(
            ch_inputs,
        )

        // Get inputs
        ch_sage_inputs = CHANNEL_GROUP_INPUTS.out.wgs_present
            .map { meta ->
                def meta_sage = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: Utils.getTumorWgsSampleName(meta),
                    normal_id: Utils.getNormalWgsSampleName(meta),
                ]
                def tumor_bam = Utils.getTumorWgsBam(meta)
                def normal_bam = Utils.getNormalWgsBam(meta)
                return [meta_sage, tumor_bam, normal_bam, "${tumor_bam}.bai", "${normal_bam}.bai"]
            }

        GERMLINE(
            ch_sage_inputs,
            ref_data_genome_fasta,
            ref_data_genome_fai,
            ref_data_genome_dict,
            ref_data_genome_version,
            ref_data_sage_known_hotspots_germline,
            ref_data_sage_actionable_panel,
            ref_data_sage_coverage_panel,
            ref_data_sage_highconf_regions,
            ref_data_ensembl_data_resources,
        )

        SOMATIC(
            ch_sage_inputs,
            ref_data_genome_fasta,
            ref_data_genome_fai,
            ref_data_genome_dict,
            ref_data_genome_version,
            ref_data_sage_known_hotspots_somatic,
            ref_data_sage_actionable_panel,
            ref_data_sage_highconf_regions,
            ref_data_ensembl_data_resources,
        )

        ch_versions = ch_versions.mix(
            GERMLINE.out.versions,
            SOMATIC.out.versions,
        )

        ch_germline_vcf_out = WorkflowOncoanalyser.restoreMeta(GERMLINE.out.vcf_filtered, ch_inputs)
        ch_germline_coverage_out = WorkflowOncoanalyser.restoreMeta(GERMLINE.out.gene_coverage, ch_inputs)
        ch_somatic_vcf_out = WorkflowOncoanalyser.restoreMeta(SOMATIC.out.vcf_filtered, ch_inputs)
        ch_somatic_tumor_bqr_out = WorkflowOncoanalyser.restoreMeta(SOMATIC.out.tumor_bqr_png, ch_inputs)
        ch_somatic_normal_bqr_out = WorkflowOncoanalyser.restoreMeta(SOMATIC.out.normal_bqr_png, ch_inputs)

    emit:
        germline_vcf       = ch_germline_vcf_out       // channel: [val(meta), sage_vcf, sage_tbi]
        germline_coverage  = ch_germline_coverage_out  // channel: [val(meta), sage_coverage]
        somatic_vcf        = ch_somatic_vcf_out        // channel: [val(meta), sage_vcf, sage_tbi]
        somatic_tumor_bqr  = ch_somatic_tumor_bqr_out  // channel: [val(meta), sage_bqr_plot]
        somatic_normal_bqr = ch_somatic_normal_bqr_out // channel: [val(meta), sage_brq_plot]

        versions           = ch_versions                    // channel: [versions.yml]
}
