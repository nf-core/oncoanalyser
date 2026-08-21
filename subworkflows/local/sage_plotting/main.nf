//
// SAGE is a precise and highly sensitive somatic SNV, MNV and small INDEL caller
//

nextflow.enable.types = true

include { SAGE_VISUALISER  } from '../../../modules/local/sage/visualiser/main'

include { FileType                    } from '../utils_nfcore_oncoanalyser_pipeline/types'
include { groupByMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { joinMeta                    } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { restoreMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { getDonorDnaSample           } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getDonorDnaSampleName       } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getDonorReduxDirAlignment   } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getDonorReduxTsvs           } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getInput                    } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getNormalDnaSample          } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getNormalDnaSampleName      } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getNormalReduxDirAlignment  } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getNormalReduxTsvs          } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSample           } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSampleName       } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorReduxDirAlignment   } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorReduxTsvs           } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { hasInput                    } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { selectCurrentOrExisting     } from '../utils_nfcore_oncoanalyser_pipeline/utils'

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
    ch_inputs_sorted = groupByMeta([
        ch_redux_dir_tumor,
        ch_redux_dir_normal,
        ch_redux_dir_donor,
        ch_purple_dir,
    ])
        .map { meta, redux_dir_tumor, redux_dir_normal, redux_dir_donor, purple_dir ->

            def redux_dir_tumor_selected = selectCurrentOrExisting(redux_dir_tumor, getInput(getTumorDnaSample(meta), FileType.REDUX_DIR))
            def redux_dir_normal_selected = selectCurrentOrExisting(redux_dir_normal, getInput(getNormalDnaSample(meta), FileType.REDUX_DIR))
            def redux_dir_donor_selected = selectCurrentOrExisting(redux_dir_donor, getInput(getDonorDnaSample(meta), FileType.REDUX_DIR))

            def (tumor_aln, tumor_idx) = getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)
            def (normal_aln, normal_idx) = getNormalReduxDirAlignment(meta, redux_dir_normal_selected)
            def (donor_aln, donor_idx) = getDonorReduxDirAlignment(meta, redux_dir_donor_selected)

            def redux_tsvs_tumor = getTumorReduxTsvs(meta, redux_dir_tumor_selected)
            def redux_tsvs_normal = getNormalReduxTsvs(meta, redux_dir_normal_selected)
            def redux_tsvs_donor = getDonorReduxTsvs(meta, redux_dir_donor_selected)
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
                selectCurrentOrExisting(purple_dir, getInput(getTumorDnaSample(meta), FileType.PURPLE_DIR)),
            ]

        }
        .branch { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx, redux_tsvs, purple_dir ->

            def has_existing = hasInput(getTumorDnaSample(meta), FileType.SAGE_PLOT_DIR)

            def tumor_dna_id = getTumorDnaSampleName(meta)
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
                key: meta.case_id,
                id: meta.case_id,
                tumor_id: getTumorDnaSampleName(meta),
            ]

            if (normal_aln) {
                meta_sage.normal_id = getNormalDnaSampleName(meta)
            }

            if (donor_aln) {
                meta_sage.donor_id = getDonorDnaSampleName(meta)
            }

            def purple_smlv_vcf = purple_dir.resolve("${getTumorDnaSampleName(meta)}.purple.somatic.vcf.gz")
            def purple_smlv_vcf_tbi = purple_dir.resolve("${getTumorDnaSampleName(meta)}.purple.somatic.vcf.gz.tbi")

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
            restoreMeta(channel.topic('sage_visualiser_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    emit:
    visualiser_dir = ch_outputs  // channel: [ meta, sage_visualiser_dir ]
}
