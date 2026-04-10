//
// PURPLE is a CNV caller that infers purity/ploidy and recovers low-confidence SVs
//

include { PURPLE } from '../../../modules/local/purple/main'

workflow PURPLE_CALLING {
    take:
    // Sample data
    ch_inputs                    // channel: [mandatory] [ meta ]
    ch_amber_dir                 // channel: [mandatory] [ meta, amber_dir ]
    ch_cobalt_dir                // channel: [mandatory] [ meta, cobalt_dir ]
    ch_esvee_dir                 // channel: [mandatory] [ meta, esvee_dir ]
    ch_pave_somatic_dir          // channel: [mandatory] [ meta, pave_dir ]
    ch_pave_germline_dir         // channel: [mandatory] [ meta, pave_dir ]
    ch_redux_tumor_dir           // channel: [optional]  [ meta, redux_dir ]

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

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources
    // channel: { meta, amber_dir, cobalt_dir, esvee_dir, pave_somatic_dir, pave_germline_dir, redux_tumor_dir }
    ch_inputs_selected = WorkflowChannels.groupByMeta(
        flatten_mode: 'singletons_only',
        ch_amber_dir,
        ch_cobalt_dir,
        ch_esvee_dir,
        ch_pave_somatic_dir,
        ch_pave_germline_dir,
        ch_redux_tumor_dir,
    )
        .map { meta, amber_dir, cobalt_dir, esvee_dir, pave_somatic_dir, pave_germline_dir, redux_tumor_dir ->

            def inputs = [:]

            inputs.meta              = meta
            inputs.amber_dir         = Inputs.preferUserProvidedInput(amber_dir, meta, Inputs.KEY.AMBER_DIR)
            inputs.cobalt_dir        = Inputs.preferUserProvidedInput(cobalt_dir, meta, Inputs.KEY.COBALT_DIR)
            inputs.esvee_dir         = Inputs.preferUserProvidedInput(esvee_dir, meta, Inputs.KEY.ESVEE_DIR)
            inputs.pave_somatic_dir  = Inputs.preferUserProvidedInput(pave_somatic_dir, meta, Inputs.KEY.PAVE_DIR_TUMOR)
            inputs.pave_germline_dir = Inputs.preferUserProvidedInput(pave_germline_dir, meta, Inputs.KEY.PAVE_DIR_NORMAL)
            inputs.redux_tumor_tsvs  = Inputs.resolveReduxTsvFiles(redux_tumor_dir, meta, samplesheet.SampleType.TUMOR)

            return inputs
        }

    // Sort inputs
    // channel: runnable: { meta, amber_dir, cobalt_dir, esvee_dir, pave_somatic_dir, pave_germline_dir, [redux_tumor_tsv, ...] }
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_inputs_selected
        .branch { inputs ->

            def has_existing = Inputs.hasExistingInput(inputs.meta, Inputs.KEY.PURPLE_DIR)

            runnable: inputs.amber_dir && inputs.cobalt_dir && !has_existing
                return inputs
            skip: true
                return inputs.meta
        }

    // Create process input channel
    // channel: [ meta_purple, amber_dir, cobalt_dir, esvee_dir, pave_somatic_dir, pave_germline_dir, [redux_tumor_tsv, ...] ]
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

            return inputs.values()
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
    )

    ch_versions = ch_versions.mix(PURPLE.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, purple_dir ]
    ch_outputs = Channel.empty()
        .mix(
            WorkflowChannels.restoreMeta(PURPLE.out.purple_dir, ch_inputs),
            channels.PlaceholderChannels.toolDir(ch_inputs_sorted.skip),
        )

    emit:
    purple_dir = ch_outputs  // channel: [ meta, purple_dir ]

    versions   = ch_versions // channel: [ versions.yml ]
}
