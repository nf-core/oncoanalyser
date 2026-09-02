//
// VIRUSBreakend and Virus Interpreter identify viral content and insertion sites
//

include { VIRUSBREAKEND    } from '../../../modules/local/virusbreakend/main'
include { VIRUSINTERPRETER } from '../../../modules/local/virusinterpreter/main'

workflow VIRUSBREAKEND_CALLING {
    take:
    // Sample data
    ch_inputs             // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor    // channel: [mandatory] [ meta, redux_dir ]
    ch_bamtools_dir_tumor // channel: [mandatory] [ meta, bamtools ]
    ch_purple             // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_fasta          // channel: [mandatory] /path/to/genome_fasta
    genome_fai            // channel: [mandatory] /path/to/genome_fai
    genome_dict           // channel: [mandatory] /path/to/genome_dict
    genome_gridss_index   // channel: [mandatory] /path/to/genome_gridss_index
    virusbreakenddb       // channel: [mandatory] /path/to/virusbreakenddb/
    virus_taxonomy_db     // channel: [mandatory] /path/to/virus_taxonomy_db
    virus_reporting_db    // channel: [mandatory] /path/to/virus_reporting_db
    virus_blocklist_db    // channel: [mandatory] /path/to/virus_blocklist_db

    // Params
    gridss_config         // channel: [optional] /path/to/gridss_config

    main:
    //
    // STEP: Handle inputs
    //
    // Select input sources then sort
    // NOTE(SW): VIRUSBreakend inputs are not allowed in the samplesheet, so aren't considered
    // channel: [ meta, tumor_aln, tumor_idx ]
    ch_inputs_sorted = ch_redux_dir_tumor
        .map { meta, redux_dir_tumor ->

            def redux_dir_tumor_selected = Utils.selectCurrentOrExisting(redux_dir_tumor, meta, Constants.INPUT.REDUX_DIR_TUMOR)
            def (tumor_aln, tumor_idx) = Utils.getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)

            return [meta, tumor_aln, tumor_idx]

        }
        .branch { meta, tumor_aln, tumor_idx ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.VIRUSINTERPRETER_DIR)
            runnable: tumor_aln && ! has_existing
            skip: true
                return meta
        }

    //
    // MODULE: VIRUSBreakend
    //
    // Create process input channel
    // channel: [ meta_virus, tumor_aln ]
    ch_virusbreakend_inputs = ch_inputs_sorted.runnable
        .map { meta, tumor_aln, _tumor_idx ->

            def meta_virus = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
            ]

            return [meta_virus, tumor_aln]
        }

    // Run process
    VIRUSBREAKEND(
        ch_virusbreakend_inputs,
        genome_fasta,
        genome_fai,
        genome_dict,
        genome_gridss_index,
        virusbreakenddb,
        gridss_config,
    )

    //
    // MODULE: Virus Interpreter
    //
    // Select input sources
    // channel: [ meta, virusbreakend_tsv, bamtools_dir_tumor, purple_dir ]
    ch_virusinterpreter_inputs_selected = WorkflowOncoanalyser.groupByMeta(
        WorkflowOncoanalyser.restoreMeta(channel.topic('virusbreakend_tsv'), ch_inputs),
        ch_bamtools_dir_tumor,
        ch_purple,
    )
        .map { meta, virusbreakend_tsv, bamtools_dir_tumor, purple_dir ->

            return [
                meta,
                virusbreakend_tsv,
                Utils.selectCurrentOrExisting(bamtools_dir_tumor, meta, Constants.INPUT.BAMTOOLS_DIR_TUMOR),
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
            ]

        }

    // Sort inputs
    // channel: [ meta, virusbreakend_tsv, bamtools_dir_tumor, purple_dir ]
    // channel: skip: [ meta ]
    ch_virusinterpreter_inputs_sorted = ch_virusinterpreter_inputs_selected
        .branch { meta, virusbreakend_tsv, bamtools_dir_tumor, purple_dir ->
            runnable: virusbreakend_tsv && bamtools_dir_tumor && purple_dir
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_virus, virusbreakend_tsv, bamtools_dir_tumor, purple_dir ]
    ch_virusinterpreter_inputs = ch_virusinterpreter_inputs_sorted.runnable
        .map { d ->

            def meta = d[0]
            def inputs = d[1..-1]

            def meta_virus = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
            ]

            return [meta_virus] + inputs
        }

    // Run process
    VIRUSINTERPRETER(
        ch_virusinterpreter_inputs,
        virus_taxonomy_db,
        virus_reporting_db,
        virus_blocklist_db,
    )

    //
    // STEP: Handle outputs
    //
    // Set outputs, restoring original meta
    // channel: [ meta, virusinterpreter_dir ]
    ch_outputs = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('virusinterpreter_dir'), ch_inputs),
            ch_virusinterpreter_inputs_sorted.skip.map { meta -> [meta, []] },
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    virusinterpreter_dir = ch_outputs // channel: [ meta, virusinterpreter_dir ]
}
