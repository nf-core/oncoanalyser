//
// SAGE is a precise and highly sensitive somatic SNV, MNV and small INDEL caller
//

import Constants
import Utils

include { SAGE_VISUALISER } from '../../../modules/local/sage/visualiser/main'

workflow SAGE_PLOTTING {
    take:
    // Sample data
    ch_inputs                   // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor          // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal         // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_donor          // channel: [mandatory] [ meta, redux_dir ]
    ch_purple_dir               // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_fasta                // channel: [mandatory] /path/to/genome_fasta
    genome_version              // channel: [mandatory] genome version
    genome_fai                  // channel: [mandatory] /path/to/genome_fai
    genome_dict                 // channel: [mandatory] /path/to/genome_dict
    sage_pon                    // channel: [mandatory] /path/to/sage_pon
    sage_known_hotspots_somatic // channel: [mandatory] /path/to/sage_known_hotspots_somatic
    sage_highconf_regions       // channel: [mandatory] /path/to/sage_highconf_regions
    ensembl_data_resources      // channel: [mandatory] /path/to/ensembl_data_resources/

    // Params
    targeted_mode               // boolean: [mandatory] Set targeted mode

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources then sort
    // channel: runnable: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, [redux_tsv, ...], purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_redux_dir_tumor,
        ch_redux_dir_normal,
        ch_redux_dir_donor,
        ch_purple_dir,
    )
        .map { meta, redux_dir_tumor, redux_dir_normal, redux_dir_donor, purple_dir ->

            def redux_dir_tumor_selected = Utils.selectCurrentOrExisting(redux_dir_tumor, meta, Constants.INPUT.REDUX_DIR_TUMOR)
            def redux_dir_normal_selected = Utils.selectCurrentOrExisting(redux_dir_normal, meta, Constants.INPUT.REDUX_DIR_NORMAL)
            def redux_dir_donor_selected = Utils.selectCurrentOrExisting(redux_dir_donor, meta, Constants.INPUT.REDUX_DIR_DONOR)

            def (tumor_bam, tumor_bai) = Utils.getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)
            def (normal_bam, normal_bai) = Utils.getNormalReduxDirAlignment(meta, redux_dir_normal_selected)
            def (donor_bam, donor_bai) = Utils.getDonorReduxDirAlignment(meta, redux_dir_donor_selected)

            def redux_tsvs_tumor = Utils.getTumorReduxTsvs(meta, redux_dir_tumor_selected)
            def redux_tsvs_normal = Utils.getNormalReduxTsvs(meta, redux_dir_normal_selected)
            def redux_tsvs_donor = Utils.getDonorReduxTsvs(meta, redux_dir_donor_selected)
            def redux_tsvs = [*redux_tsvs_tumor, *redux_tsvs_normal, *redux_tsvs_donor].findAll{ it.exists() }

            return [
                meta,
                tumor_bam,
                tumor_bai,
                normal_bam,
                normal_bai,
                donor_bam,
                donor_bai,
                redux_tsvs,
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
            ]

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
            ]

            if (normal_bam) {
                meta_sage.normal_id = Utils.getNormalDnaSampleName(meta)
            }

            if (donor_bam) {
                meta_sage.donor_id = Utils.getDonorDnaSampleName(meta)
            }

            def purple_smlv_vcf = purple_dir.resolve("${Utils.getTumorDnaSampleName(meta)}.purple.somatic.vcf.gz")
            def purple_smlv_vcf_tbi = purple_dir.resolve("${Utils.getTumorDnaSampleName(meta)}.purple.somatic.vcf.gz.tbi")

            return [
                meta_sage,
                tumor_bam,
                normal_bam,
                donor_bam,
                tumor_bai,
                normal_bai,
                donor_bai,
                redux_tsvs,
                purple_smlv_vcf,
                purple_smlv_vcf_tbi,
            ]

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
        targeted_mode,
    )

    ch_versions = ch_versions.mix(SAGE_VISUALISER.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, sage_visualiser_dir ]
    ch_outputs = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(SAGE_VISUALISER.out.sage_visualiser_dir, ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    visualiser_dir = ch_outputs  // channel: [ meta, sage_plot_dir ]

    versions       = ch_versions // channel: [ versions.yml ]
}
