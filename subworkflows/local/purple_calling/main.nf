//
// PURPLE is a CNV caller that infers purity/ploidy and recovers low-confidence SVs
//

import Constants
import Inputs

include { PURPLE } from '../../../modules/local/purple/main'

workflow PURPLE_CALLING {
    take:
    // Sample data
    ch_inputs                    // channel: [mandatory] [ meta ]
    ch_amber                     // channel: [mandatory] [ meta, amber_dir ]
    ch_cobalt                    // channel: [mandatory] [ meta, cobalt_dir ]
    ch_pave_somatic              // channel: [mandatory] [ meta, pave_dir ]
    ch_pave_germline             // channel: [mandatory] [ meta, pave_dir ]
    ch_esvee                     // channel: [mandatory] [ meta, esvee_dir ]

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
    germline_amp_del_freq        // channel: [optional]  /path/to/germline_amp_del_freq
    target_region_bed            // channel: [optional]  /path/to/target_region_bed
    target_region_ratios         // channel: [optional]  /path/to/target_region_ratios
    target_region_msi_indels     // channel: [optional]  /path/to/target_region_msi_indels

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources
    // channel: [ meta, amber_dir, cobalt_dir, pave_somatic_dir, pave_germline_dir, esvee_dir ]
    ch_inputs_selected = WorkflowOncoanalyser.groupByMeta(
        ch_amber,
        ch_cobalt,
        ch_pave_somatic,
        ch_pave_germline,
        ch_esvee,
    )
        .map { d ->

            def meta = d[0]

            // NOTE(SW): avoiding further complexity with loops etc

            def inputs = [
                Inputs.overrideWithExistingInput(d[1], meta, Constants.INPUT.AMBER_DIR),
                Inputs.overrideWithExistingInput(d[2], meta, Constants.INPUT.COBALT_DIR),
                Inputs.overrideWithExistingInput(d[3], meta, Constants.INPUT.PAVE_DIR_TUMOR),
                Inputs.overrideWithExistingInput(d[4], meta, Constants.INPUT.PAVE_DIR_NORMAL),
                Inputs.overrideWithExistingInput(d[5], meta, Constants.INPUT.ESVEE_DIR),
            ]

            return [meta, *inputs]
        }

    // Sort inputs
    // channel: runnable: [ meta, amber_dir, cobalt_dir, pave_somatic_dir, pave_germline_dir, esvee_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_inputs_selected
        .branch { d ->
            def meta = d[0]
            def amber_dir = d[1]
            def cobalt_dir = d[2]

            def has_existing = Inputs.hasExistingInput(meta, Constants.INPUT.PURPLE_DIR)

            runnable: amber_dir && cobalt_dir && !has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_purple, amber_dir, cobalt_dir, pave_somatic_dir, pave_germline_dir, esvee_dir ]
    ch_purple_inputs = ch_inputs_sorted.runnable
        .map { d ->

            def meta = d[0]
            def inputs = d[1..-1]

            def meta_purple = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Inputs.getTumorDnaSampleName(meta),
            ]

            if (Inputs.hasNormalDna(meta)) {
                meta_purple.normal_id = Inputs.getNormalDnaSampleName(meta)
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
        germline_amp_del_freq,
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
            PlaceholderChannels.toolDir(ch_inputs_sorted.skip),
        )

    emit:
    purple_dir = ch_outputs  // channel: [ meta, purple_dir ]

    versions   = ch_versions // channel: [ versions.yml ]
}
