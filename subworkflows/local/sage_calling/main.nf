//
// SAGE is a precise and highly sensitive somatic SNV, MNV and small INDEL caller
//

import Constants
import Utils

import java.nio.channels.Channel

include { SAGE_GERMLINE } from '../../../modules/local/sage/germline/main'
include { SAGE_SOMATIC  } from '../../../modules/local/sage/somatic/main'

workflow SAGE_CALLING {
    take:
    // Sample data
    ch_inputs                    // channel: [mandatory] [ meta ]
    ch_tumor_bam                 // channel: [mandatory] [ meta, bam, bai ]
    ch_normal_bam                // channel: [mandatory] [ meta, bam, bai ]
    ch_donor_bam                 // channel: [mandatory] [ meta, bam, bai ]
    ch_tumor_tsv                 // channel: [mandatory] [ meta, dup_freq_tsv, jitter_tsv, ms_tsv ]
    ch_normal_tsv                // channel: [mandatory] [ meta, dup_freq_tsv, jitter_tsv, ms_tsv ]
    ch_donor_tsv                 // channel: [mandatory] [ meta, dup_freq_tsv, jitter_tsv, ms_tsv ]

    // Reference data
    genome_fasta                 // channel: [mandatory] /path/to/genome_fasta
    genome_version               // channel: [mandatory] genome version
    genome_fai                   // channel: [mandatory] /path/to/genome_fai
    genome_dict                  // channel: [mandatory] /path/to/genome_dict
    sage_known_hotspots_somatic  // channel: [mandatory] /path/to/sage_known_hotspots_somatic
    sage_known_hotspots_germline // channel: [optional]  /path/to/sage_known_hotspots_germline
    sage_actionable_panel        // channel: [mandatory] /path/to/sage_actionable_panel
    sage_coverage_panel          // channel: [mandatory] /path/to/sage_coverage_panel
    sage_highconf_regions        // channel: [mandatory] /path/to/sage_highconf_regions
    segment_mappability          // channel: [mandatory] /path/to/segment_mappability
    driver_gene_panel            // channel: [mandatory] /path/to/driver_gene_panel
    ensembl_data_resources       // channel: [mandatory] /path/to/ensembl_data_resources/

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Sort inputs
    // channel: runnable: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, [redux_tsv, ...] ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_tumor_bam,
        ch_normal_bam,
        ch_donor_bam,
        ch_tumor_tsv,
        ch_normal_tsv,
        ch_donor_tsv,
    )
        .map { meta,
            tumor_bam,  tumor_bai,
            normal_bam, normal_bai,
            donor_bam,  donor_bai,

            tumor_dup_freq_tsv,  tumor_jitter_tsv,  tumor_ms_tsv,
            normal_dup_freq_tsv, normal_jitter_tsv, normal_ms_tsv,
            donor_dup_freq_tsv,  donor_jitter_tsv,  donor_ms_tsv ->

            def redux_tsv_list = [
                tumor_dup_freq_tsv  ?: Utils.getInput(meta, Constants.INPUT.REDUX_DUP_FREQ_TSV_TUMOR),
                tumor_jitter_tsv    ?: Utils.getInput(meta, Constants.INPUT.REDUX_JITTER_TSV_TUMOR),
                tumor_ms_tsv        ?: Utils.getInput(meta, Constants.INPUT.REDUX_MS_TSV_TUMOR),

                normal_dup_freq_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_DUP_FREQ_TSV_NORMAL),
                normal_jitter_tsv   ?: Utils.getInput(meta, Constants.INPUT.REDUX_JITTER_TSV_NORMAL),
                normal_ms_tsv       ?: Utils.getInput(meta, Constants.INPUT.REDUX_MS_TSV_NORMAL),

                donor_dup_freq_tsv  ?: Utils.getInput(meta, Constants.INPUT.REDUX_DUP_FREQ_TSV_DONOR),
                donor_jitter_tsv    ?: Utils.getInput(meta, Constants.INPUT.REDUX_JITTER_TSV_DONOR),
                donor_ms_tsv        ?: Utils.getInput(meta, Constants.INPUT.REDUX_MS_TSV_DONOR),
            ]

            redux_tsv_list = redux_tsv_list.findAll{ it != [] }

            return [
                meta,

                Utils.selectCurrentOrExisting(tumor_bam, meta, Constants.INPUT.BAM_REDUX_DNA_TUMOR),
                tumor_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_TUMOR),

                Utils.selectCurrentOrExisting(normal_bam, meta, Constants.INPUT.BAM_REDUX_DNA_NORMAL),
                normal_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_NORMAL),

                Utils.selectCurrentOrExisting(donor_bam, meta, Constants.INPUT.BAM_REDUX_DNA_DONOR),
                donor_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_DONOR),

                redux_tsv_list,
            ]
        }
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, redux_tsvs ->
            runnable: tumor_bam
            skip: true
                return meta
        }

    //
    // MODULE: SAGE germline
    //
    // Select inputs that are eligible to run
    // channel: runnable: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, [redux_tsv, ...] ]
    // channel: skip: [ meta ]
    ch_inputs_germline_sorted = ch_inputs_sorted.runnable
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, redux_tsvs ->
            def has_tumor_normal = tumor_bam && normal_bam
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SAGE_VCF_NORMAL)

            runnable: has_tumor_normal && !has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_sage, tumor_bam, normal_bam, tumor_bai, normal_bai, [redux_tsv, ...] ]
    ch_sage_germline_inputs = ch_inputs_germline_sorted.runnable
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, redux_tsvs ->

            def meta_sage = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Utils.getTumorDnaSampleName(meta),
                normal_id: Utils.getNormalDnaSampleName(meta),
            ]

            return [meta_sage, tumor_bam, normal_bam, tumor_bai, normal_bai, redux_tsvs]
        }

    // Run process
    SAGE_GERMLINE(
        ch_sage_germline_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        sage_known_hotspots_germline,
        sage_actionable_panel,
        sage_coverage_panel,
        sage_highconf_regions,
        ensembl_data_resources,
    )

    ch_versions = ch_versions.mix(SAGE_GERMLINE.out.versions)

    //
    // MODULE: SAGE somatic
    //
    // Select inputs that are eligible to run
    // channel: runnable: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, [redux_tsv, ...] ]
    // channel: skip: [ meta ]
    ch_inputs_somatic_sorted = ch_inputs_sorted.runnable
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, redux_tsvs ->
            def has_tumor = tumor_bam
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SAGE_VCF_TUMOR)

            runnable: has_tumor && !has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: tumor/normal: [ meta_sage, tumor_bam, normal_bam, donor_bam, tumor_bai, normal_bai, donor_bai, [redux_tsv, ...] ]
    // channel: tumor only: [ meta_sage, tumor_bam, [], tumor_bai, [], [redux_tsv, ...] ]
    ch_sage_somatic_inputs = ch_inputs_somatic_sorted.runnable
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, redux_tsvs ->

            def meta_sage = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Utils.getTumorDnaSampleName(meta),
                donor_id: Utils.getDonorDnaSampleName(meta),
            ]

            if (normal_bam) {
                meta_sage.normal_id = Utils.getNormalDnaSampleName(meta)
            }

            if (donor_bam) {
                meta_sage.donor_id = Utils.getDonorDnaSampleName(meta)
            }

            return [meta_sage, tumor_bam, normal_bam, donor_bam, tumor_bai, normal_bai, donor_bai, redux_tsvs]
        }

    // Run process
    SAGE_SOMATIC(
        ch_sage_somatic_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        sage_known_hotspots_somatic,
        sage_actionable_panel,
        sage_coverage_panel,
        sage_highconf_regions,
        ensembl_data_resources,
    )

    ch_versions = ch_versions.mix(SAGE_SOMATIC.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, sage_vcf, sage_tbi ]
    ch_somatic_vcf_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(SAGE_SOMATIC.out.vcf, ch_inputs),
            ch_inputs_somatic_sorted.skip.map { meta -> [meta, [], []] },
            ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
        )

    // channel: [ meta, sage_vcf, sage_tbi ]
    ch_germline_vcf_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(SAGE_GERMLINE.out.vcf, ch_inputs),
            ch_inputs_germline_sorted.skip.map { meta -> [meta, [], []] },
            ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
        )

    // channel: [ meta, sage_dir ]
    ch_somatic_dir = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(SAGE_SOMATIC.out.sage_dir, ch_inputs),
            ch_inputs_somatic_sorted.skip.map { meta -> [meta, []] },
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    // channel: [ meta, sage_dir ]
    ch_germline_dir = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(SAGE_GERMLINE.out.sage_dir, ch_inputs),
            ch_inputs_germline_sorted.skip.map { meta -> [meta, []] },
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    germline_vcf = ch_germline_vcf_out // channel: [ meta, sage_vcf, sage_tbi ]
    somatic_vcf  = ch_somatic_vcf_out  // channel: [ meta, sage_vcf, sage_tbi ]
    germline_dir = ch_germline_dir     // channel: [ meta, sage_dir ]
    somatic_dir  = ch_somatic_dir      // channel: [ meta, sage_dir ]

    versions     = ch_versions         // channel: [ versions.yml ]
}
