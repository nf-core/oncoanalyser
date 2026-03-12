//
// PURPLE is a CNV caller that infers purity/ploidy and recovers low-confidence SVs
//

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
    // channel: { meta, amber_dir, cobalt_dir, pave_somatic_dir, pave_germline_dir, esvee_dir }
    ch_inputs_selected = WorkflowOncoanalyser.groupByMeta(
        ch_amber,
        ch_cobalt,
        ch_pave_somatic,
        ch_pave_germline,
        ch_esvee,
    )
        .map { meta, amber_dir, cobalt_dir, pave_somatic_dir, pave_germline_dir, esvee_dir ->

            def inputs = [:]

            inputs.meta              = meta
            inputs.amber_dir         = Inputs.preferUserProvidedInput(amber_dir, meta, Constants.INPUT.AMBER_DIR)
            inputs.cobalt_dir        = Inputs.preferUserProvidedInput(cobalt_dir, meta, Constants.INPUT.COBALT_DIR)
            inputs.pave_somatic_dir  = Inputs.preferUserProvidedInput(pave_somatic_dir, meta, Constants.INPUT.PAVE_DIR_TUMOR)
            inputs.pave_germline_dir = Inputs.preferUserProvidedInput(pave_germline_dir, meta, Constants.INPUT.PAVE_DIR_NORMAL)
            inputs.esvee_dir         = Inputs.preferUserProvidedInput(esvee_dir, meta, Constants.INPUT.ESVEE_DIR)

            return inputs
        }

    // Sort inputs
    // channel: runnable: { meta, amber_dir, cobalt_dir, pave_somatic_dir, pave_germline_dir, esvee_dir }
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_inputs_selected
        .branch { inputs ->

            def has_existing = Inputs.hasExistingInput(inputs.meta, Constants.INPUT.PURPLE_DIR)

            runnable: inputs.amber_dir && inputs.cobalt_dir && !has_existing
                return inputs
            skip: true
                return inputs.meta
        }

    // Create process input channel
    // channel: [ meta_purple, amber_dir, cobalt_dir, pave_somatic_dir, pave_germline_dir, esvee_dir ]
    ch_purple_inputs = ch_inputs_sorted.runnable
        .map { inputs ->

            def meta = inputs.meta

            def meta_purple = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Inputs.getTumorDnaSampleName(meta),
            ]

            if (Inputs.hasNormalDna(meta)) {
                meta_purple.normal_id = Inputs.getNormalDnaSampleName(meta)
            }

            inputs.meta = meta_purple

            return inputs
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
