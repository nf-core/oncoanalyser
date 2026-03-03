//
// SAGE is a precise and highly sensitive somatic SNV, MNV and small INDEL caller
//

import java.nio.channels.Channel

include { SAGE_VISUALISER } from '../../../modules/local/sage/visualiser/main'

workflow SAGE_PLOTTING {
    take:
    // Sample data
    ch_inputs                    // channel: [mandatory] [ meta ]
    ch_tumor_bam                 // channel: [mandatory] [ meta, bam, bai ]
    ch_normal_bam                // channel: [mandatory] [ meta, bam, bai ]
    ch_donor_bam                 // channel: [mandatory] [ meta, bam, bai ]
    ch_tumor_tsv                 // channel: [mandatory] [ meta, redux_tsv, ... ]
    ch_normal_tsv                // channel: [mandatory] [ meta, redux_tsv, ... ]
    ch_donor_tsv                 // channel: [mandatory] [ meta, redux_tsv, ... ]
    ch_purple                    // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_fasta                 // channel: [mandatory] /path/to/genome_fasta
    genome_version               // channel: [mandatory] genome version
    genome_fai                   // channel: [mandatory] /path/to/genome_fai
    genome_dict                  // channel: [mandatory] /path/to/genome_dict
    sage_pon                     // channel: [mandatory] /path/to/sage_pon
    sage_known_hotspots_somatic  // channel: [mandatory] /path/to/sage_known_hotspots_somatic
    sage_highconf_regions        // channel: [mandatory] /path/to/sage_highconf_regions
    ensembl_data_resources       // channel: [mandatory] /path/to/ensembl_data_resources/

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Sort inputs
    // channel: runnable: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, [redux_tsv, ...], purple_dir ]
    // channel: skip: [ meta ]

    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_tumor_bam, ch_tumor_tsv,
        ch_normal_bam, ch_normal_tsv,
        ch_donor_bam, ch_donor_tsv,
        ch_purple,
    )
        .map { meta,
            tumor_bam , tumor_bai , tumor_bqr_tsv , tumor_dup_freq_tsv , tumor_jitter_tsv , tumor_ms_tsv,
            normal_bam, normal_bai, normal_bqr_tsv, normal_dup_freq_tsv, normal_jitter_tsv, normal_ms_tsv,
            donor_bam , donor_bai , donor_bqr_tsv , donor_dup_freq_tsv , donor_jitter_tsv , donor_ms_tsv,
            purple_dir ->

            tumor_bam = Utils.overrideWithExistingInput(tumor_bam, meta, Constants.INPUT.BAM_REDUX_DNA_TUMOR)
            tumor_bai = tumor_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_TUMOR)

            normal_bam = Utils.overrideWithExistingInput(normal_bam, meta, Constants.INPUT.BAM_REDUX_DNA_NORMAL)
            normal_bai = normal_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_NORMAL)

            donor_bam = Utils.overrideWithExistingInput(donor_bam, meta, Constants.INPUT.BAM_REDUX_DNA_DONOR)
            donor_bai = donor_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_DONOR)

            def redux_tsvs = [
                tumor_bqr_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_BQR_TSV_TUMOR),
                tumor_jitter_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_JITTER_TSV_TUMOR),
                tumor_ms_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_MS_TSV_TUMOR),

                normal_bqr_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_BQR_TSV_NORMAL),
                normal_jitter_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_JITTER_TSV_NORMAL),
                normal_ms_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_MS_TSV_NORMAL),

                donor_bqr_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_BQR_TSV_DONOR),
                donor_jitter_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_JITTER_TSV_DONOR),
                donor_ms_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_MS_TSV_DONOR),
            ]

            redux_tsvs = redux_tsvs.findAll { it != [] }

            purple_dir = Utils.overrideWithExistingInput(purple_dir, meta, Constants.INPUT.PURPLE_DIR)

            return [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, redux_tsvs, purple_dir ]
        }
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, redux_tsvs, purple_dir ->
            runnable: tumor_bam && purple_dir
            skip: true
            return meta
        }

    // Create process input channel
    // channel: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, [redux_tsv, ...], purple_smlv_vcf ]
    ch_sage_plotting_inputs = ch_inputs_sorted.runnable
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, redux_tsvs, purple_dir ->

            def meta_sage = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Utils.getTumorDnaSampleName(meta),
                normal_id: normal_bam ? Utils.getNormalDnaSampleName(meta) : null,
                donor_id: donor_bam ? Utils.getDonorDnaSampleName(meta) : null,
            ]

            def purple_smlv_vcf = Utils.getPurpleSomaticVcf(meta, purple_dir)
            def purple_smlv_vcf_tbi = Utils.getPurpleSomaticVcfTbi(meta, purple_dir)

            return [meta_sage, tumor_bam, normal_bam, donor_bam, tumor_bai, normal_bai, donor_bai, redux_tsvs, purple_smlv_vcf, purple_smlv_vcf_tbi]
        }

    // Run process
    SAGE_VISUALISER(
        ch_sage_plotting_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        sage_pon,
        sage_known_hotspots_somatic,
        sage_highconf_regions,
        ensembl_data_resources,
    )

    ch_versions = ch_versions.mix(SAGE_VISUALISER.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, sage_dir ]
    ch_visualiser_dir_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(SAGE_VISUALISER.out.sage_vis_dir, ch_inputs),
            PlaceholderChannels.toolDir(ch_inputs_sorted.skip),
        )

    emit:
    visualiser_dir  = ch_visualiser_dir_out // channel: [ meta, sage_dir ]
    versions        = ch_versions           // channel: [ versions.yml ]
}
