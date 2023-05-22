//
// SAGE is a precise and highly sensitive somatic SNV, MNV and small INDEL caller.
//

include { SAGE_GERMLINE as GERMLINE } from '../../modules/local/sage/germline/main'
include { SAGE_SOMATIC as SOMATIC   } from '../../modules/local/sage/somatic/main'

workflow SAGE_CALLING {
    take:
        // Sample data
        ch_inputs                             // channel: [val(meta)]

        // Reference data
        ref_data_genome_fasta                 //    file: /path/to/genome_fasta
        ref_data_genome_fai                   //    file: /path/to/genome_fai
        ref_data_genome_dict                  //    file: /path/to/genome_dict
        ref_data_genome_version               //     val: genome version
        ref_data_sage_known_hotspots_germline //    file: /path/to/sage_known_hotspots_germline
        ref_data_sage_known_hotspots_somatic  //    file: /path/to/sage_known_hotspots_somatic
        ref_data_sage_actionable_panel        //    file: /path/to/sage_actionable_panel
        ref_data_sage_coverage_panel          //    file: /path/to/sage_coverage_panel
        ref_data_sage_highconf_regions        //    file: /path/to/sage_highconf_regions
        ref_data_segment_mappability          //    file: /path/to/segment_mappability
        ref_data_driver_gene_panel            //    file: /path/to/driver_gene_panel
        ref_data_ensembl_data_resources       //    file: /path/to/ensembl_data_resources/
        ref_data_ensembl_data_resources       //    file: /path/to/ensembl_data_resources/

        // Params
        run_config

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Get inputs
        ch_sage_inputs = ch_inputs
            .map { meta ->
                def meta_sage = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: Utils.getTumorSampleName(meta, run_config.mode),
                ]

                def tumor_bam = Utils.getTumorBam(meta, run_config.mode)

                def normal_bam = []
                def normal_bai = []

                if (run_config.type == Constants.RunType.TUMOR_NORMAL) {

                    assert [Constants.RunMode.WGS, Constants.RunMode.WGTS].contains(run_config.mode)

                    meta_sage.normal_id = Utils.getNormalWgsSampleName(meta)
                    normal_bam = Utils.getNormalWgsBam(meta)
                    normal_bai = "${normal_bam}.bai"

                }

                return [meta_sage, tumor_bam, normal_bam, "${tumor_bam}.bai", normal_bai]
            }

        //
        // MODULE: SAGE germline
        //
        ch_germline_vcf_out = Channel.empty()
        ch_germline_coverage_out = Channel.empty()
        if (run_config.type == Constants.RunType.TUMOR_NORMAL) {

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

            ch_versions = ch_versions.mix(GERMLINE.out.versions)
            ch_germline_vcf_out = WorkflowOncoanalyser.restoreMeta(GERMLINE.out.vcf_filtered, ch_inputs)
            ch_germline_coverage_out = WorkflowOncoanalyser.restoreMeta(GERMLINE.out.gene_coverage, ch_inputs)
        }

        //
        // MODULE: SAGE somatic
        //
        SOMATIC(
            ch_sage_inputs,
            ref_data_genome_fasta,
            ref_data_genome_fai,
            ref_data_genome_dict,
            ref_data_genome_version,
            ref_data_sage_known_hotspots_somatic,
            ref_data_sage_actionable_panel,
            ref_data_sage_coverage_panel,
            ref_data_sage_highconf_regions,
            ref_data_ensembl_data_resources,
        )

        ch_versions = ch_versions.mix(SOMATIC.out.versions)
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
