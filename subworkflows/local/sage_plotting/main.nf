//
// SAGE is a precise and highly sensitive somatic SNV, MNV and small INDEL caller
//

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
    // Select input sources then sort
    // channel: runnable: [ meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx, [redux_tsv, ...], purple_dir ]
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

            def (tumor_aln, tumor_idx) = Utils.getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)
            def (normal_aln, normal_idx) = Utils.getNormalReduxDirAlignment(meta, redux_dir_normal_selected)
            def (donor_aln, donor_idx) = Utils.getDonorReduxDirAlignment(meta, redux_dir_donor_selected)

            def redux_tsvs_tumor = Utils.getTumorReduxTsvs(meta, redux_dir_tumor_selected)
            def redux_tsvs_normal = Utils.getNormalReduxTsvs(meta, redux_dir_normal_selected)
            def redux_tsvs_donor = Utils.getDonorReduxTsvs(meta, redux_dir_donor_selected)
            def redux_tsvs = (redux_tsvs_tumor + redux_tsvs_normal + redux_tsvs_donor).findAll{ tsv -> tsv.exists() }

            return [
                meta,
                tumor_aln,
                tumor_idx,
                normal_aln,
                normal_idx,
                donor_aln,
                donor_idx,
                redux_tsvs,
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
            ]

        }
        .branch { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx, redux_tsvs, purple_dir ->

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SAGE_PLOT_DIR_TUMOR)

            def tumor_dna_id = Utils.getTumorDnaSampleName(meta)
            def has_smlv_vcf = purple_dir ? purple_dir.resolve("${tumor_dna_id}.purple.somatic.vcf.gz").exists() : false

            runnable: tumor_aln && has_smlv_vcf && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx, [redux_tsv, ...], purple_smlv_vcf ]
    ch_sage_plotting_inputs = ch_inputs_sorted.runnable
        .map { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx, redux_tsvs, purple_dir ->

            def meta_sage = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Utils.getTumorDnaSampleName(meta),
            ]

            if (normal_aln) {
                meta_sage.normal_id = Utils.getNormalDnaSampleName(meta)
            }

            if (donor_aln) {
                meta_sage.donor_id = Utils.getDonorDnaSampleName(meta)
            }

            def purple_smlv_vcf = purple_dir.resolve("${Utils.getTumorDnaSampleName(meta)}.purple.somatic.vcf.gz")
            def purple_smlv_vcf_tbi = purple_dir.resolve("${Utils.getTumorDnaSampleName(meta)}.purple.somatic.vcf.gz.tbi")

            return [
                meta_sage,
                tumor_aln,
                normal_aln,
                donor_aln,
                tumor_idx,
                normal_idx,
                donor_idx,
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

    // Set outputs, restoring original meta
    // channel: [ meta, sage_visualiser_dir ]
    ch_outputs = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('sage_visualiser_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    visualiser_dir = ch_outputs  // channel: [ meta, sage_visualiser_dir ]
}
