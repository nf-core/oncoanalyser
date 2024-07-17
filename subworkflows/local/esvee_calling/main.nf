//
// ESVEE detects structural variants, and reports breakends and breakpoints.
//

import Constants
import Utils

import java.nio.channels.Channel

include { ESVEE_PREP            } from '../../../modules/local/esvee/svprep/main'
include { ESVEE_ASSEMBLE        } from '../../../modules/local/esvee/assemble/main'
include { ESVEE_DEPTH_ANNOTATOR } from '../../../modules/local/esvee/depth_annotator/main'
include { ESVEE_CALL            } from '../../../modules/local/esvee/call/main'

workflow ESVEE_CALLING {
    take:
    // Sample data
    ch_inputs                // channel: [mandatory] [ meta ]
    ch_tumor_bam             // channel: [mandatory] [ meta, bam, bai ]
    ch_normal_bam            // channel: [mandatory] [ meta, bam, bai ]

    // Reference data
    genome_fasta             // channel: [mandatory] /path/to/genome_fasta
    genome_version           // channel: [mandatory] genome version
    genome_fai               // channel: [mandatory] /path/to/genome_fai
    genome_dict              // channel: [mandatory] /path/to/genome_dict
    genome_gridss_index      // channel: [mandatory] /path/to/genome_gridss_index
    gridss_blocklist         // channel: [mandatory] /path/to/gridss_blocklist
    sv_prep_blocklist        // channel: [mandatory] /path/to/sv_prep_blocklist
    known_fusions            // channel: [mandatory] /path/to/known_fusions
    decoy_sequences          // channel: [mandatory] /path/to/deocy_sequences // TODO: make optional
    pon_breakends            // channel: [mandatory] /path/to/pon_sgl
    pon_breakpoints          // channel: [mandatory] /path/to/pon_sv
    repeatmasker_annotations // channel: [mandatory] /path/to/repeatmasker_annotations


    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources and sort
    // channel: runnable_tn: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai ]
    // channel: runnable_to: [ meta, tumor_bam, tumor_bai ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_tumor_bam,
        ch_normal_bam,
    )
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            return [
                meta,
                Utils.selectCurrentOrExisting(tumor_bam, meta, Constants.INPUT.BAM_MARKDUPS_DNA_TUMOR),
                tumor_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_TUMOR),
                Utils.selectCurrentOrExisting(normal_bam, meta, Constants.INPUT.BAM_MARKDUPS_DNA_NORMAL),
                normal_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_NORMAL),
            ]
        }
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.ESVEE_VCF_TUMOR)

            runnable_tn: tumor_bam && normal_bam && !has_existing
            runnable_to: tumor_bam && !has_existing
                return [meta, tumor_bam, tumor_bai]
            skip: true
                return meta
        }

    //
    // MODULE: Esvee prep
    //
    // Create process input channel
    // channel: [ meta_svprep, bam_tumor, bai_tumor, [] ]
    ch_sv_prep_inputs = Channel.empty()
        .mix(
            ch_inputs_sorted.runnable_to.map { [*it, [], []] },
            ch_inputs_sorted.runnable_tn,
        )
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->

            def meta_esvee_prep = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Utils.getTumorDnaSampleName(meta),
                normal_id: Utils.getNormalDnaSampleName(meta),
            ]

            return [meta_esvee_prep, tumor_bam, normal_bam, tumor_bai, normal_bai]
        }

    // Run process
    ESVEE_PREP(
        ch_sv_prep_inputs,
        genome_fasta,
        genome_version,
        sv_prep_blocklist,
        known_fusions,
    )

    ch_versions = ch_versions.mix(ESVEE_PREP.out.versions)


    //
    // MODULE: ESVEE assemble reads
    //
    // Create process input channel
    // channel: [ meta_esvee, sv_prep_dir ]
    ch_assemble_inputs =
        WorkflowOncoanalyser.groupByMeta(
            ch_sv_prep_inputs,
            ESVEE_PREP.out.tumor_prep_bam,
            ESVEE_PREP.out.normal_prep_bam,
            ESVEE_PREP.out.junctions_tsv,
            ESVEE_PREP.out.fragment_lengths_tsv,
        )
        .map {
            meta, tumor_bam, tumor_bai, normal_bam, normal_bai,
            tumor_prep_bam, normal_prep_bam, junctions_tsv, fragment_lengths_tsv ->
            return [meta, tumor_prep_bam, normal_prep_bam, junctions_tsv, fragment_lengths_tsv]
        }

    // Run process
    ESVEE_ASSEMBLE(
        ch_assemble_inputs,
        genome_fasta,
        genome_fai,
        genome_dict,
        genome_version,
        decoy_sequences,
    )

    ch_versions = ch_versions.mix(ESVEE_ASSEMBLE.out.versions)


    // MODULE: ESVEE annotated reference sample depth
    //
    // Create process input channel
    // channel: [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, assemble_dir]
    ch_depth_annotator_inputs =
        WorkflowOncoanalyser.groupByMeta(
            ch_sv_prep_inputs,
            ESVEE_ASSEMBLE.out.raw_vcf,
        )
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, raw_vcf ->
            return [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, raw_vcf]
        }

    // Run process
    ESVEE_DEPTH_ANNOTATOR(
        ch_depth_annotator_inputs,
        genome_fasta,
        genome_version,
    )

    ch_versions = ch_versions.mix(ESVEE_DEPTH_ANNOTATOR.out.versions)


    //
    // MODULE: ESVEE call somatic structural variants
    //
    // Create process input channel
    // channel: [meta, ref_depth_vcf, fragment_lengths_tsv]
    ch_call_inputs =
        WorkflowOncoanalyser.groupByMeta(
            ch_sv_prep_inputs,
            ESVEE_DEPTH_ANNOTATOR.out.ref_depth_vcf,
            ESVEE_PREP.out.fragment_lengths_tsv,
        )
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, ref_depth_vcf, fragment_lengths_tsv ->
            return [meta, ref_depth_vcf, fragment_lengths_tsv]
        }

    ESVEE_CALL(
        ch_call_inputs,
        genome_fasta,
        genome_version,
        pon_breakends,
        pon_breakpoints,
        known_fusions,
        repeatmasker_annotations,
    )

    ch_versions = ch_versions.mix(ESVEE_CALL.out.versions)


    // Set outputs, restoring original meta
    // channel: [ meta, esvee_vcf ]
    ch_outputs = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ESVEE_CALL.out.somatic_vcf, ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    vcf      = ch_outputs  // channel: [ meta, vcf ]

    versions = ch_versions // channel: [ versions.yml ]
}
