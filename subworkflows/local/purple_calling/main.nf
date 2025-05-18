//
// PURPLE is a CNV caller that infers purity/ploidy and recovers low-confidence SVs
//

import Constants
import Utils

include { PURPLE } from '../../../modules/local/purple/main'

workflow PURPLE_CALLING {
    take:
    // Sample data
    ch_inputs                    // channel: [mandatory] [ meta ]
    ch_amber                     // channel: [mandatory] [ meta, amber_dir ]
    ch_cobalt                    // channel: [mandatory] [ meta, cobalt_dir ]
    ch_smlv_somatic              // channel: [mandatory] [ meta, pave_vcf ]
    ch_smlv_germline             // channel: [mandatory] [ meta, pave_vcf ]
    ch_sv_somatic                // channel: [mandatory] [ meta, esvee_vcf, esvee_tbi ]
    ch_sv_germline               // channel: [mandatory] [ meta, esvee_vcf, esvee_tbi ]

    // Reference data
    genome_fasta                 // channel: [mandatory] /path/to/genome_fasta
    genome_version               // channel: [mandatory] genome version
    genome_fai                   // channel: [mandatory] /path/to/genome_fai
    genome_dict                  // channel: [mandatory] /path/to/genome_dict
    gc_profile                   // channel: [mandatory] /path/to/gc_profile
    sage_known_hotspots_somatic  // channel: [mandatory] /path/to/sage_known_hotspots_somatic
    sage_known_hotspots_germline // channel: [optional]  /path/to/sage_known_hotspots_germline
    driver_gene_panel            // channel: [mandatory] /path/to/driver_gene_panel
    ensembl_data_resources       // channel: [mandatory] /path/to/ensembl_data_resources/
    purple_germline_del          // channel: [optional]  /path/to/purple_germline_del
    target_region_bed            // channel: [optional]  /path/to/target_region_bed
    target_region_ratios         // channel: [optional]  /path/to/target_region_ratios
    target_region_msi_indels     // channel: [optional]  /path/to/target_region_msi_indels

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources
    // channel: [ meta, amber_dir, cobalt_dir, sv_somatic_vcf, sv_somatic_tbi, sv_germline_vcf, sv_germline_tbi, smlv_somatic_vcf, smlv_germline_vcf ]
    ch_inputs_selected = WorkflowOncoanalyser.groupByMeta(
        ch_amber,
        ch_cobalt,
        ch_sv_somatic,
        ch_sv_germline,
        ch_smlv_somatic,
        ch_smlv_germline,
    )
        .map { d ->

            def meta = d[0]

            // NOTE(SW): avoiding further complexity with loops etc

            def inputs = [
                Utils.selectCurrentOrExisting(d[1], meta, Constants.INPUT.AMBER_DIR),
                Utils.selectCurrentOrExisting(d[2], meta, Constants.INPUT.COBALT_DIR),
                Utils.selectCurrentOrExisting(d[3], meta, Constants.INPUT.ESVEE_VCF_TUMOR),
                Utils.selectCurrentOrExisting(d[4], meta, Constants.INPUT.ESVEE_VCF_TUMOR_TBI),
                Utils.selectCurrentOrExisting(d[5], meta, Constants.INPUT.ESVEE_VCF_NORMAL),
                Utils.selectCurrentOrExisting(d[6], meta, Constants.INPUT.ESVEE_VCF_NORMAL_TBI),
                Utils.selectCurrentOrExisting(d[7], meta, Constants.INPUT.PAVE_VCF_TUMOR),
                Utils.selectCurrentOrExisting(d[8], meta, Constants.INPUT.PAVE_VCF_NORMAL),
            ]

            return [meta, *inputs]
        }

    // Sort inputs
    // channel: runnable: [ meta, amber_dir, cobalt_dir, sv_somatic_vcf, sv_somatic_tbi, sv_germline_vcf, sv_germline_tbi, smlv_somatic_vcf, smlv_germline_vcf ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_inputs_selected
        .branch { d ->
            def meta = d[0]
            def amber_dir = d[1]
            def cobalt_dir = d[2]

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.PURPLE_DIR)

            runnable: amber_dir && cobalt_dir && !has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_purple, amber_dir, cobalt_dir, sv_somatic_vcf, sv_somatic_tbi, sv_germline_vcf, sv_germline_tbi, smlv_somatic_vcf, smlv_germline_vcf ]
    ch_purple_inputs = ch_inputs_sorted.runnable
        .map { d ->

            def meta = d[0]
            def inputs = d[1..-1]

            def meta_purple = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Utils.getTumorDnaSampleName(meta),
            ]

            if (Utils.hasNormalDna(meta)) {
                meta_purple.normal_id = Utils.getNormalDnaSampleName(meta)
            }

            return [meta_purple, *inputs]

        }

    // Run process
    PURPLE(
        ch_purple_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        gc_profile,
        sage_known_hotspots_somatic,
        sage_known_hotspots_germline,
        driver_gene_panel,
        ensembl_data_resources,
        purple_germline_del,
        target_region_bed,
        target_region_ratios,
        target_region_msi_indels,
    )

    ch_versions = ch_versions.mix(PURPLE.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, purple_dir ]
    ch_outputs = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(PURPLE.out.purple_dir, ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    purple_dir = ch_outputs  // channel: [ meta, purple_dir ]

    versions   = ch_versions // channel: [ versions.yml ]
}
