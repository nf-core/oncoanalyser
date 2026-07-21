//
// ESVEE detects structural variants, and reports breakends and breakpoints.
//

import Constants
import Utils

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
    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Select input sources then sort
    // channel: runnable: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_redux_dir_tumor,
        ch_redux_dir_normal,
    )
        .map { meta, redux_dir_tumor, redux_dir_normal ->

            def redux_dir_tumor_selected = Utils.selectCurrentOrExisting(redux_dir_tumor, meta, Constants.INPUT.REDUX_DIR_TUMOR)
            def redux_dir_normal_selected = Utils.selectCurrentOrExisting(redux_dir_normal, meta, Constants.INPUT.REDUX_DIR_NORMAL)

            def (tumor_bam, tumor_bai) = Utils.getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)
            def (normal_bam, normal_bai) = Utils.getNormalReduxDirAlignment(meta, redux_dir_normal_selected)

            return [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
        }
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.ESVEE_DIR)

            runnable_tn: tumor_bam && normal_bam && ! has_existing
            runnable_to: tumor_bam && ! has_existing
                return [meta, tumor_bam, tumor_bai]
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_esvee, tumor_bam, tumor_bai, normal_bam, normal_bai ]
    ch_esvee_inputs = Channel.empty()
        .mix(
            ch_inputs_sorted.runnable_tn,
            ch_inputs_sorted.runnable_to.map { [*it, [], []] },
        )
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->

            def meta_esvee = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Utils.getTumorDnaSampleName(meta),
            ]

            if (normal_bam) {
                meta_esvee.normal_id = Utils.getNormalDnaSampleName(meta)
            }

            return [meta_esvee, tumor_bam, tumor_bai, normal_bam, normal_bai]
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

    ch_versions = ch_versions.mix(ESVEE.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, esvee_dir ]
    ch_outputs = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ESVEE.out.esvee_dir, ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    esvee_dir = ch_outputs  // channel: [ meta, esvee_dir ]

    versions  = ch_versions // channel: [ versions.yml ]
}
