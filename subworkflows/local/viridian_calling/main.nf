//
// Viridian calls oncovirus presence and integration sites
//

include { VIRIDIAN } from '../../../modules/local/viridian/main'

workflow VIRIDIAN_CALLING {
    take:
    // Sample data
    ch_inputs             // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor    // channel: [mandatory] [ meta, redux_dir ]
    ch_esvee_dir          // channel: [mandatory] [ meta, esvee_dir ]
    ch_bamtools_dir_tumor // channel: [mandatory] [ meta, bamtools_dir ]
    ch_purple_dir         // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_fasta          // channel: [mandatory] /path/to/genome_fasta
    genome_fai            // channel: [mandatory] /path/to/genome_fai
    genome_dict           // channel: [mandatory] /path/to/genome_dict
    viridian_resources    // channel: [mandatory] /path/to/viridian_resources/

    main:
    //
    // STEP: Handle inputs
    //
    // Select input sources then sort
    // channel: runnable: [ meta, tumor_aln, tumor_idx, esvee_dir, bamtools_dir_tumor, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_redux_dir_tumor,
        ch_esvee_dir,
        ch_bamtools_dir_tumor,
        ch_purple_dir,
    )
        .map { meta, redux_dir_tumor, esvee_dir, bamtools_dir_tumor, purple_dir ->

            def redux_dir_tumor_selected = Utils.selectCurrentOrExisting(redux_dir_tumor, meta, Constants.INPUT.REDUX_DIR_TUMOR)
            def (tumor_aln, tumor_idx) = Utils.getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)

            return [
                meta,
                tumor_aln,
                tumor_idx,
                Utils.selectCurrentOrExisting(esvee_dir, meta, Constants.INPUT.ESVEE_DIR),
                Utils.selectCurrentOrExisting(bamtools_dir_tumor, meta, Constants.INPUT.BAMTOOLS_DIR_TUMOR),
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
            ]

        }
        .branch { meta, tumor_aln, tumor_idx, esvee_dir, bamtools_dir_tumor, purple_dir ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.VIRIDIAN_DIR)
            runnable: tumor_aln && esvee_dir && bamtools_dir_tumor && purple_dir && ! has_existing
            skip: true
                return meta
        }

    //
    // MODULE: Viridian
    //
    // Create process input channel
    // channel: [ meta_viridian, tumor_aln, tumor_idx, esvee_dir, bamtools_dir_tumor, purple_dir ]
    ch_viridian_inputs = ch_inputs_sorted.runnable
        .map { d ->

            def meta = d[0]
            def inputs = d[1..-1]

            def meta_viridian = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
            ]

            return [meta_viridian] + inputs
        }

    // Run process
    VIRIDIAN(
        ch_viridian_inputs,
        genome_fasta,
        genome_fai,
        genome_dict,
        viridian_resources,
    )

    //
    // STEP: Handle outputs
    //
    // Set outputs, restoring original meta
    // channel: [ meta, viridian_dir ]
    ch_outputs = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('viridian_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    viridian_dir = ch_outputs // channel: [ meta, viridian_dir ]
}
