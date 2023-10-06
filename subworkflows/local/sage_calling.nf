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

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Sort inputs
        // channel: [ meta ]
        ch_inputs_sorted = ch_inputs.branch { meta ->
            runnable: Utils.hasTumorDnaBam(meta)
            skip: true
        }

        // Create general process input channel
        // channel: tumor/normal: [ meta_sage, tbam, nbam, tbai, nbai ]
        // channel: tumor only:   [ meta_sage, tbam, [], tbai, [] ]
        ch_sage_inputs = ch_inputs_sorted.runnable
            .map { meta ->

                // NOTE(SW): germline only is not currently supported
                assert Utils.hasTumorDnaBam(meta)

                def meta_sage = [
                    // Both are required. Reminder:
                    //   * key: channel element grouping
                    //   * id: task tag
                    key: meta.group_id,
                    id: meta.group_id,
                    tumor_id: Utils.getTumorDnaSampleName(meta)
                ]

                def data = []

                if (Utils.hasNormalDnaBam(meta)) {

                    meta_sage.normal_id = Utils.getNormalDnaSampleName(meta)
                    meta_sage.sample_type = 'tumor_normal'

                    data = [
                        meta_sage,
                        Utils.getTumorDnaBam(meta),
                        Utils.getNormalDnaBam(meta),
                        Utils.getTumorDnaBai(meta),
                        Utils.getNormalDnaBai(meta),
                    ]

                } else if (Utils.hasTumorDnaBam(meta)) {

                    meta_sage.sample_type = 'tumor_only'

                    data = [
                        meta_sage,
                        Utils.getTumorDnaBam(meta),
                        [],
                        Utils.getTumorDnaBai(meta),
                        [],
                    ]

                } else {
                    assert false
                }

                return data
            }

        //
        // MODULE: SAGE germline
        //
        // Create germline input channel
        // channel: [ meta_sage, tbam, nbam, tbai, nbai ]
        ch_sage_germline_inputs = ch_sage_inputs
            .filter {
                def meta_sage = it[0]
                return meta_sage.sample_type == 'tumor_normal'
            }

        // Run process
        GERMLINE(
            ch_sage_germline_inputs,
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

        //
        // MODULE: SAGE somatic
        //
        // Run process
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

        // Set outputs, restoring original meta
        // NOTE(SW): look to reduce repetitiveness here
        // channel: [ meta, sage_vcf, sage_tbi ]
        ch_somatic_vcf_out = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(SOMATIC.out.vcf_filtered, ch_inputs),
                ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
            )

        // channel: [ meta, sage_vcf, sage_tbi ]
        ch_germline_vcf_out = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(GERMLINE.out.vcf_filtered, ch_inputs),
                ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
            )

        // channel: [ meta, sage_dir ]
        ch_somatic_dir = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(SOMATIC.out.sage_dir, ch_inputs),
                ch_inputs_sorted.skip.map { meta -> [meta, []] },
            )

        // channel: [ meta, sage_dir ]
        ch_germline_dir = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(GERMLINE.out.sage_dir, ch_inputs),
                ch_inputs_sorted.skip.map { meta -> [meta, []] },
            )

    emit:
        germline_vcf = ch_germline_vcf_out // channel: [ meta, sage_vcf, sage_tbi ]
        somatic_vcf  = ch_somatic_vcf_out  // channel: [ meta, sage_vcf, sage_tbi ]
        germline_dir = ch_germline_dir     // channel: [ meta, sage_dir ]
        somatic_dir  = ch_somatic_dir      // channel: [ meta, sage_dir ]

        versions           = ch_versions   // channel: [ versions.yml ]
}
