//
// ESVEE detects structural variants, and reports breakends and breakpoints.
//

include { ESVEE } from '../../../modules/local/esvee/main'

workflow ESVEE_CALLING {
    take:

    // Sample data
    ch_inputs                // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor       // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal      // channel: [mandatory] [ meta, redux_dir ]

    // Reference data
    genome_fasta             // channel: [mandatory] /path/to/genome_fasta
    genome_version           // channel: [mandatory] genome version
    genome_fai               // channel: [mandatory] /path/to/genome_fai
    genome_dict              // channel: [mandatory] /path/to/genome_dict
    genome_img               // channel: [optional]  /path/to/genome_img
    known_fusions            // channel: [mandatory] /path/to/known_fusions
    pon_breakends            // channel: [mandatory] /path/to/pon_sgl
    pon_breakpoints          // channel: [mandatory] /path/to/pon_sv
    decoy_sequences_image    // channel: [mandatory] /path/to/decoy_sequences_image
    repeatmasker_annotations // channel: [mandatory] /path/to/repeatmasker_annotations
    unmap_regions            // channel: [mandatory] /path/to/unmap_regions
    target_regions_bed       // channel: [optional]  /path/to/target_regions_bed

    // Params
    sequencing_platform      // string:  [mandatory] sequencing platform

    main:
    // Select input sources then sort
    // channel: runnable: [ meta, tumor_aln, tumor_idx, normal_aln, normal_idx ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_redux_dir_tumor,
        ch_redux_dir_normal,
    )
        .map { meta, redux_dir_tumor, redux_dir_normal ->

            def redux_dir_tumor_selected = Utils.selectCurrentOrExisting(redux_dir_tumor, meta, Constants.INPUT.REDUX_DIR_TUMOR)
            def redux_dir_normal_selected = Utils.selectCurrentOrExisting(redux_dir_normal, meta, Constants.INPUT.REDUX_DIR_NORMAL)

            def (tumor_aln, tumor_idx) = Utils.getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)
            def (normal_aln, normal_idx) = Utils.getNormalReduxDirAlignment(meta, redux_dir_normal_selected)

            return [meta, tumor_aln, tumor_idx, normal_aln, normal_idx]
        }
        .branch { meta, tumor_aln, tumor_idx, normal_aln, normal_idx ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.ESVEE_DIR)

            runnable_tn: tumor_aln && normal_aln && ! has_existing
            runnable_to: tumor_aln && ! has_existing
                return [meta, tumor_aln, tumor_idx]
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_esvee, tumor_aln, tumor_idx, normal_aln, normal_idx ]
    ch_esvee_inputs = channel.empty()
        .mix(
            ch_inputs_sorted.runnable_tn,
            ch_inputs_sorted.runnable_to.map { d -> d + [[], []] },
        )
        .map { meta, tumor_aln, tumor_idx, normal_aln, normal_idx ->

            def meta_esvee = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Utils.getTumorDnaSampleName(meta),
            ]

            if (normal_aln) {
                meta_esvee.normal_id = Utils.getNormalDnaSampleName(meta)
            }

            return [meta_esvee, tumor_aln, tumor_idx, normal_aln, normal_idx]
        }

    // Run process
    ESVEE(
        ch_esvee_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        genome_img,
        pon_breakends,
        pon_breakpoints,
        decoy_sequences_image,
        known_fusions,
        repeatmasker_annotations,
        unmap_regions,
        target_regions_bed,
        sequencing_platform,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, esvee_dir ]
    ch_outputs = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('esvee_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    esvee_dir = ch_outputs // channel: [ meta, esvee_dir ]
}
