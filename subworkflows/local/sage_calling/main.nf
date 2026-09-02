//
// SAGE is a precise and highly sensitive somatic SNV, MNV and small INDEL caller
//

nextflow.enable.types = true

include { SAGE_GERMLINE } from '../../../modules/local/sage/germline/main'
include { SAGE_SOMATIC  } from '../../../modules/local/sage/somatic/main'

include { getDonorReduxDirAlignments  } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { getDonorReduxTsvs           } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { getNormalReduxDirAlignment  } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { getNormalReduxTsvs          } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { getTumorReduxDirAlignment   } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { getTumorReduxTsvs           } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { getDonorDnaSample           } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getDonorDnaSampleNames      } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getInput                    } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalDnaSample          } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalDnaSampleName      } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSample           } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSampleName       } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasInput                    } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { FileType                    } from '../utils_nfcore_oncoanalyser_pipeline/types_enums'
include { groupByMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { joinMeta                    } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { restoreMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { selectCurrentOrExisting     } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow SAGE_CALLING {
    take:
    // Sample data
    ch_inputs                   : Channel<Map>              // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor          : Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal         : Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_donor          : Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, redux_dir ]

    // Reference data
    genome_fasta                : Channel<Path>             // channel: [mandatory] /path/to/genome_fasta
    genome_version              : Channel<String>           // channel: [mandatory] genome version
    genome_fai                  : Channel<Path>             // channel: [mandatory] /path/to/genome_fai
    genome_dict                 : Channel<Path>             // channel: [mandatory] /path/to/genome_dict
    sage_pon                    : Channel<Path>             // channel: [mandatory] /path/to/sage_pon
    sage_known_hotspots_somatic : Channel<Path>             // channel: [mandatory] /path/to/sage_known_hotspots_somatic
    sage_known_hotspots_germline: Channel<Path>?            // channel: [optional]  /path/to/sage_known_hotspots_germline
    sage_highconf_regions       : Channel<Path>             // channel: [mandatory] /path/to/sage_highconf_regions
    segment_mappability         : Channel<Path>             // channel: [mandatory] /path/to/segment_mappability
    driver_gene_panel           : Channel<Path>             // channel: [mandatory] /path/to/driver_gene_panel
    ensembl_data_resources      : Channel<Path>             // channel: [mandatory] /path/to/ensembl_data_resources/
    gnomad_resource             : Channel<Path>             // channel: [mandatory] /path/to/gnomad_resource

    // Params
    sequencing_platform         : String                    // string:  [mandatory] sequencing platform
    targeted_mode               : Boolean                   // boolean: [mandatory] Set targeted mode
    enable_germline             : Boolean                   // boolean: [mandatory] Enable germline mode

    main:
    //
    // STEP: Handle inputs
    //
    // Select input sources then sort
    // channel: [ meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx, [redux_tsv, ...] ]
    ch_inputs_sorted = groupByMeta([
        ch_redux_dir_tumor,
        ch_redux_dir_normal,
        ch_redux_dir_donor,
    ])
        .map { meta, redux_dir_tumor, redux_dir_normal, redux_dir_donor ->

            def redux_dir_tumor_selected = selectCurrentOrExisting(redux_dir_tumor, getInput(getTumorDnaSample(meta), FileType.REDUX_DIR))
            def redux_dir_normal_selected = selectCurrentOrExisting(redux_dir_normal, getInput(getNormalDnaSample(meta), FileType.REDUX_DIR))
            def redux_dir_donor_selected = selectCurrentOrExisting(redux_dir_donor, getInput(getDonorDnaSample(meta), FileType.REDUX_DIR))

            def (tumor_aln, tumor_idx) = getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)
            def (normal_aln, normal_idx) = getNormalReduxDirAlignment(meta, redux_dir_normal_selected)
            def donor_alignments = getDonorReduxDirAlignments(meta, redux_dir_donor_selected)
            def donor_alns = donor_alignments.collect { aln, idx -> aln }
            def donor_idxs = donor_alignments.collect { aln, idx -> idx }

            def redux_tsvs_tumor = getTumorReduxTsvs(meta, redux_dir_tumor_selected)
            def redux_tsvs_normal = getNormalReduxTsvs(meta, redux_dir_normal_selected)
            def redux_tsvs_donor = getDonorReduxTsvs(meta, redux_dir_donor_selected)
            def redux_tsvs = (redux_tsvs_tumor + redux_tsvs_normal + redux_tsvs_donor).findAll{ tsv -> tsv.exists() }

            return [meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_alns, donor_idxs, redux_tsvs]
        }
        .branch { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_alns, donor_idxs, redux_tsvs ->
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
        .branch { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_alns, donor_idxs, redux_tsvs ->
            def has_tumor_normal = tumor_aln && normal_aln
            def has_existing = hasInput(getNormalDnaSample(meta), FileType.SAGE_DIR)

            runnable: has_tumor_normal && ! has_existing && enable_germline
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_sage, tumor_aln, tumor_idx, normal_aln, normal_idx, [redux_tsv, ...] ]
    ch_sage_germline_inputs = ch_inputs_germline_sorted.runnable
        .map { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, _donor_alns, _donor_idxs, redux_tsvs ->

            def meta_sage = record(
                key: meta.case_id,
                id: meta.case_id,
                tumor_id: getTumorDnaSampleName(meta),
                normal_id: getNormalDnaSampleName(meta),
            )

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
        .branch { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_alns, donor_idxs, redux_tsvs ->
            def has_tumor = tumor_aln
            def has_existing = hasInput(getTumorDnaSample(meta), FileType.SAGE_DIR)

            runnable: has_tumor && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: tumor/normal: [ meta_sage, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx, [redux_tsv, ...] ]
    // channel: tumor only: [ meta_sage, tumor_aln, tumor_idx, [], [], [], [], [redux_tsv, ...] ]
    ch_sage_somatic_inputs = ch_inputs_somatic_sorted.runnable
        .map { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_alns, donor_idxs, redux_tsvs ->

            def meta_sage = record(
                key: meta.case_id,
                id: meta.case_id,
                tumor_id: getTumorDnaSampleName(meta),
                normal_id: normal_aln ? getNormalDnaSampleName(meta) : null,
                donor_ids: donor_alns ? getDonorDnaSampleNames(meta) : null,
            )

            return [meta_sage, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_alns, donor_idxs, redux_tsvs]
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
            restoreMeta(channel.topic('sage_somatic_dir'), ch_inputs),
            ch_inputs_somatic_sorted.skip.map { meta -> [meta, null] },
            ch_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    // channel: [ meta, sage_dir ]
    ch_outputs_germline = channel.empty()
        .mix(
            restoreMeta(channel.topic('sage_germline_dir'), ch_inputs),
            ch_inputs_germline_sorted.skip.map { meta -> [meta, null] },
            ch_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    emit:
    germline_dir = ch_outputs_germline // channel: [ meta, sage_dir ]
    somatic_dir  = ch_outputs_somatic  // channel: [ meta, sage_dir ]
}
