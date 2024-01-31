//
// SAGE is a precise and highly sensitive somatic SNV, MNV and small INDEL caller
//

import Constants
import Utils

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
        ch_inputs_sorted = ch_inputs
            .branch { meta ->
                runnable: Utils.hasTumorDnaBam(meta)
                skip: true
            }

        //
        // MODULE: SAGE germline
        //
        // Select inputs that are eligible to run
        // channel: [ meta ]
        ch_inputs_germline_sorted = ch_inputs_sorted.runnable
            .branch { meta ->
                def has_tumor_normal = Utils.hasTumorDnaBam(meta) && Utils.hasNormalDnaBam(meta)
                def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SAGE_VCF_NORMAL)

                runnable: has_tumor_normal && !has_existing
                skip: true
            }

        // Create process input channel
        // channel: [ meta_sage, tbam, nbam, tbai, nbai ]
        ch_sage_germline_inputs = ch_inputs_germline_sorted.runnable
            .map { meta ->

                def meta_sage = [
                    key: meta.group_id,
                    id: meta.group_id,
                    tumor_id: Utils.getTumorDnaSampleName(meta),
                    normal_id: Utils.getNormalDnaSampleName(meta),
                ]

                data = [
                    meta_sage,
                    Utils.getTumorDnaBam(meta),
                    Utils.getNormalDnaBam(meta),
                    Utils.getTumorDnaBai(meta),
                    Utils.getNormalDnaBai(meta),
                ]

                return data

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
        // Select inputs that are eligible to run
        // channel: [ meta ]
        ch_inputs_somatic_sorted = ch_inputs_sorted.runnable
            .branch { meta ->
                def has_tumor = Utils.hasTumorDnaBam(meta)
                def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SAGE_VCF_TUMOR)

                runnable: has_tumor && !has_existing
                skip: true
            }

        // Create process input channel
        // channel: tumor/normal: [ meta_sage, tbam, nbam, tbai, nbai ]
        // channel: tumor only: [ meta_sage, tbam, [], tbai, [] ]
        ch_sage_somatic_inputs = ch_inputs_somatic_sorted.runnable
            .map { meta ->

                def meta_sage = [
                    key: meta.group_id,
                    id: meta.group_id,
                    tumor_id: Utils.getTumorDnaSampleName(meta),
                ]

                def data = []
                if (Utils.hasNormalDnaBam(meta)) {

                    meta_sage.normal_id = Utils.getNormalDnaSampleName(meta)

                    data = [
                        meta_sage,
                        Utils.getTumorDnaBam(meta),
                        Utils.getNormalDnaBam(meta),
                        Utils.getTumorDnaBai(meta),
                        Utils.getNormalDnaBai(meta),
                    ]

                } else {

                    data = [
                        meta_sage,
                        Utils.getTumorDnaBam(meta),
                        [],
                        Utils.getTumorDnaBai(meta),
                        [],
                    ]

                }

                return data

            }

        // Run process
        SOMATIC(
            ch_sage_somatic_inputs,
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
        // channel: [ meta, sage_vcf, sage_tbi ]
        ch_somatic_vcf_out = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(SOMATIC.out.vcf, ch_inputs),
                ch_inputs_somatic_sorted.skip.map { meta -> [meta, [], []] },
                ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
            )

        // channel: [ meta, sage_vcf, sage_tbi ]
        ch_germline_vcf_out = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(GERMLINE.out.vcf, ch_inputs),
                ch_inputs_germline_sorted.skip.map { meta -> [meta, [], []] },
                ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
            )

        // channel: [ meta, sage_dir ]
        ch_somatic_dir = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(SOMATIC.out.sage_dir, ch_inputs),
                ch_inputs_somatic_sorted.skip.map { meta -> [meta, []] },
                ch_inputs_sorted.skip.map { meta -> [meta, []] },
            )

        // channel: [ meta, sage_dir ]
        ch_germline_dir = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(GERMLINE.out.sage_dir, ch_inputs),
                ch_inputs_germline_sorted.skip.map { meta -> [meta, []] },
                ch_inputs_sorted.skip.map { meta -> [meta, []] },
            )

    emit:
        germline_vcf = ch_germline_vcf_out // channel: [ meta, sage_vcf, sage_tbi ]
        somatic_vcf  = ch_somatic_vcf_out  // channel: [ meta, sage_vcf, sage_tbi ]
        germline_dir = ch_germline_dir     // channel: [ meta, sage_dir ]
        somatic_dir  = ch_somatic_dir      // channel: [ meta, sage_dir ]

        versions     = ch_versions         // channel: [ versions.yml ]
}
