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
    ch_redux_dir_tumor           // channel: [optional]  [ meta, redux_dir ]

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
    target_regions_bed           // channel: [optional]  /path/to/target_regions_bed

    main:
    // Select input sources then sort
    // channel: runnable: [ meta, amber_dir, cobalt_dir, esvee_dir, pave_somatic_dir, pave_germline_dir, [redux_tsv_tumor, ...] ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_amber_dir,
        ch_cobalt_dir,
        ch_esvee_dir,
        ch_pave_somatic_dir,
        ch_pave_germline_dir,
        ch_redux_dir_tumor,
    )
        .map { meta, amber_dir, cobalt_dir, esvee_dir, pave_somatic_dir, pave_germline_dir, redux_dir_tumor ->

            def tumor_dir_selected = Utils.selectCurrentOrExisting(redux_dir_tumor, meta, Constants.INPUT.REDUX_DIR_TUMOR)
            def tumor_tsvs = Utils.getTumorReduxTsvs(meta, tumor_dir_selected)

            return [
                meta,
                Utils.selectCurrentOrExisting(amber_dir, meta, Constants.INPUT.AMBER_DIR),
                Utils.selectCurrentOrExisting(cobalt_dir, meta, Constants.INPUT.COBALT_DIR),
                Utils.selectCurrentOrExisting(esvee_dir, meta, Constants.INPUT.ESVEE_DIR),
                Utils.selectCurrentOrExisting(pave_somatic_dir, meta, Constants.INPUT.PAVE_DIR_TUMOR),
                Utils.selectCurrentOrExisting(pave_germline_dir, meta, Constants.INPUT.PAVE_DIR_NORMAL),
                tumor_tsvs.findAll { f -> f.exists() },
            ]
        }
        .branch { meta, amber_dir, cobalt_dir, esvee_dir, pave_somatic_dir, pave_germline_dir, redux_tsvs_tumor ->

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.PURPLE_DIR)

            runnable: amber_dir && cobalt_dir && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_purple, amber_dir, cobalt_dir, esvee_dir, pave_somatic_dir, pave_germline_dir, [redux_tsv_tumor, ...] ]
    ch_purple_inputs = ch_inputs_sorted.runnable
        .map { meta, amber_dir, cobalt_dir, esvee_dir, pave_somatic_dir, pave_germline_dir, redux_tsvs_tumor ->

            def meta_purple = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Utils.getTumorDnaSampleName(meta),
            ]

            if (Utils.hasNormalDna(meta)) {
                meta_purple.normal_id = Utils.getNormalDnaSampleName(meta)
            }

            return [meta_purple, amber_dir, cobalt_dir, esvee_dir, pave_somatic_dir, pave_germline_dir, redux_tsvs_tumor]
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
        target_regions_bed,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, purple_dir ]
    ch_outputs = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('purple_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    purple_dir = ch_outputs // channel: [ meta, purple_dir ]
}
