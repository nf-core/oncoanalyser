//
// SAGE is a precise and highly sensitive somatic SNV, MNV and small INDEL caller
//

include { SAGE_GERMLINE as GERMLINE } from '../../modules/local/sage/germline/main'
include { SAGE_SOMATIC as SOMATIC   } from '../../modules/local/sage/somatic/main'

workflow SAGE_CALLING {
    take:
        // Sample data
        ch_inputs                    // channel: [mandatory] [ meta ]

        // Reference data
        genome_fasta                 // channel: [mandatory] /path/to/genome_fasta
        genome_version               // channel: [mandatory] genome version
        genome_fai                   // channel: [mandatory] /path/to/genome_fai
        genome_dict                  // channel: [mandatory] /path/to/genome_dict
        sage_known_hotspots_germline // channel: [optional]  /path/to/sage_known_hotspots_germline
        sage_known_hotspots_somatic  // channel: [mandatory] /path/to/sage_known_hotspots_somatic
        sage_actionable_panel        // channel: [mandatory] /path/to/sage_actionable_panel
        sage_coverage_panel          // channel: [mandatory] /path/to/sage_coverage_panel
        sage_highconf_regions        // channel: [mandatory] /path/to/sage_highconf_regions
        segment_mappability          // channel: [mandatory] /path/to/segment_mappability
        driver_gene_panel            // channel: [mandatory] /path/to/driver_gene_panel
        ensembl_data_resources       // channel: [mandatory] /path/to/ensembl_data_resources/
        ensembl_data_resources       // channel: [mandatory] /path/to/ensembl_data_resources/

        // Params
        run_config                   // channel: [mandatory] run configuration

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Get inputs
        // channel: [ meta_sage, tumor_bam, normal_bam, tumor_bai, normal_bai ]
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
        // channel: [ meta, sage_vcf, sage_tbi ]
        ch_germline_vcf_out = Channel.empty()
        // channel: [ meta, sage_coverage ]
        ch_germline_coverage_out = Channel.empty()
        if (run_config.type == Constants.RunType.TUMOR_NORMAL) {

            GERMLINE(
                ch_sage_inputs,
                genome_fasta,
                genome_version,
                genome_fai,
                genome_dict,
                sage_known_hotspots_germline,
                sage_actionable_panel,
                sage_coverage_panel,
                sage_highconf_regions,
                ensembl_data_resources,
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
            genome_fasta,
            genome_version,
            genome_fai,
            genome_dict,
            sage_known_hotspots_somatic,
            sage_actionable_panel,
            sage_coverage_panel,
            sage_highconf_regions,
            ensembl_data_resources,
        )

        ch_versions = ch_versions.mix(SOMATIC.out.versions)
        ch_somatic_vcf_out = WorkflowOncoanalyser.restoreMeta(SOMATIC.out.vcf_filtered, ch_inputs)
        ch_somatic_tumor_bqr_out = WorkflowOncoanalyser.restoreMeta(SOMATIC.out.tumor_bqr_png, ch_inputs)
        ch_somatic_normal_bqr_out = WorkflowOncoanalyser.restoreMeta(SOMATIC.out.normal_bqr_png, ch_inputs)

    emit:
        germline_vcf       = ch_germline_vcf_out       // channel: [ meta, sage_vcf, sage_tbi ]
        germline_coverage  = ch_germline_coverage_out  // channel: [ meta, sage_coverage ]
        somatic_vcf        = ch_somatic_vcf_out        // channel: [ meta, sage_vcf, sage_tbi ]
        somatic_tumor_bqr  = ch_somatic_tumor_bqr_out  // channel: [ meta, sage_bqr_plot ]
        somatic_normal_bqr = ch_somatic_normal_bqr_out // channel: [ meta, sage_brq_plot ]

        versions           = ch_versions               // channel: [ versions.yml ]
}
