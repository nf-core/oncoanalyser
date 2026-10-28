//
// SAGE is a precise and highly sensitive somatic SNV, MNV and small INDEL caller
//

include { SAGE_GERMLINE } from '../../../modules/local/sage/germline/main'
include { SAGE_SOMATIC  } from '../../../modules/local/sage/somatic/main'

workflow SAGE_CALLING {
    take:
    // Sample data
    ch_inputs                    // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor           // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal          // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_donor           // channel: [mandatory] [ meta, redux_dir ]

    // Reference data
    genome_fasta                 // channel: [mandatory] /path/to/genome_fasta
    genome_version               // channel: [mandatory] genome version
    genome_fai                   // channel: [mandatory] /path/to/genome_fai
    genome_dict                  // channel: [mandatory] /path/to/genome_dict
    sage_pon                     // channel: [mandatory] /path/to/sage_pon
    sage_known_hotspots_somatic  // channel: [mandatory] /path/to/sage_known_hotspots_somatic
    sage_known_hotspots_germline // channel: [optional]  /path/to/sage_known_hotspots_germline
    sage_highconf_regions        // channel: [mandatory] /path/to/sage_highconf_regions
    segment_mappability          // channel: [mandatory] /path/to/segment_mappability
    driver_gene_panel            // channel: [mandatory] /path/to/driver_gene_panel
    ensembl_data_resources       // channel: [mandatory] /path/to/ensembl_data_resources/
    gnomad_resource              // channel: [mandatory] /path/to/gnomad_resource

    // Params
    sequencing_platform          // string:  [mandatory] sequencing platform
    targeted_mode                // boolean: [mandatory] Set targeted mode
    enable_germline              // boolean: [mandatory] Enable germline mode

    main:
    //
    // STEP: Handle inputs
    //
    // Select input sources then sort
    // channel: [ meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx, [redux_tsv, ...] ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_redux_dir_tumor,
        ch_redux_dir_normal,
        ch_redux_dir_donor,
    )
        .map { meta, redux_dir_tumor, redux_dir_normal, redux_dir_donor ->

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

            return [meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx, redux_tsvs]
        }
        .branch { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx, redux_tsvs ->
            runnable: tumor_aln
            skip: true
                return meta
        }

    //
    // MODULE: SAGE germline
    //
    // Select inputs that are eligible to run
    // channel: runnable: [ meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx, [redux_tsv, ...] ]
    // channel: skip: [ meta ]
    ch_inputs_germline_sorted = ch_inputs_sorted.runnable
        .branch { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx, redux_tsvs ->
            def has_tumor_normal = tumor_aln && normal_aln
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SAGE_DIR_NORMAL)

            runnable: has_tumor_normal && ! has_existing && enable_germline
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_sage, tumor_aln, tumor_idx, normal_aln, normal_idx, [redux_tsv, ...] ]
    ch_sage_germline_inputs = ch_inputs_germline_sorted.runnable
        .map { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, _donor_aln, _donor_idx, redux_tsvs ->

            def meta_sage = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Utils.getTumorDnaSampleName(meta),
                normal_id: Utils.getNormalDnaSampleName(meta),
            ]

            return [meta_sage, tumor_aln, tumor_idx, normal_aln, normal_idx, redux_tsvs]
        }

    // Run process
    SAGE_GERMLINE(
        ch_sage_germline_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        sage_known_hotspots_germline,
        sage_highconf_regions,
        driver_gene_panel,
        ensembl_data_resources,
        sequencing_platform,
        targeted_mode,
    )

    //
    // MODULE: SAGE somatic
    //
    // Select inputs that are eligible to run
    // channel: runnable: { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx, [redux_tsv, ...] }
    // channel: skip: [ meta ]
    ch_inputs_somatic_sorted = ch_inputs_sorted.runnable
        .branch { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx, redux_tsvs ->
            def has_tumor = tumor_aln
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.SAGE_DIR_TUMOR)

            runnable: has_tumor && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: tumor/normal: [ meta_sage, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx, [redux_tsv, ...] ]
    // channel: tumor only: [ meta_sage, tumor_aln, tumor_idx, [], [], [], [], [redux_tsv, ...] ]
    ch_sage_somatic_inputs = ch_inputs_somatic_sorted.runnable
        .map { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx, redux_tsvs ->

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

            return [meta_sage, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx, redux_tsvs]
        }

    // Run process
    SAGE_SOMATIC(
        ch_sage_somatic_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        sage_pon,
        sage_known_hotspots_somatic,
        sage_highconf_regions,
        driver_gene_panel,
        ensembl_data_resources,
        gnomad_resource,
        sequencing_platform,
        targeted_mode,
    )

    //
    // STEP: Handle outputs
    //
    // Set outputs, restoring original meta
    // channel: [ meta, sage_dir ]
    ch_outputs_somatic = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('sage_somatic_dir'), ch_inputs),
            ch_inputs_somatic_sorted.skip.map { meta -> [meta, []] },
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    // channel: [ meta, sage_dir ]
    ch_outputs_germline = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('sage_germline_dir'), ch_inputs),
            ch_inputs_germline_sorted.skip.map { meta -> [meta, []] },
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    germline_dir = ch_outputs_germline // channel: [ meta, sage_dir ]
    somatic_dir  = ch_outputs_somatic  // channel: [ meta, sage_dir ]
}
