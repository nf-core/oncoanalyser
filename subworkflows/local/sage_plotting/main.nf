//
// SAGE is a precise and highly sensitive somatic SNV, MNV and small INDEL caller
//

include { SAGE_VISUALISER } from '../../../modules/local/sage/visualiser/main'

include { groupByMeta                } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { joinMeta                   } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { restoreMeta                } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { getDonorDnaReduxDir        } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getDonorDnaSampleName      } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getDonorReduxDirAlignment  } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getDonorReduxTsvs          } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getNormalDnaReduxDir       } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getNormalDnaSampleName     } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getNormalReduxDirAlignment } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getNormalReduxTsvs         } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getPurpleDir               } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getTumorDnaReduxDir        } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getTumorDnaSampleName      } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getTumorReduxDirAlignment  } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getTumorReduxTsvs          } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { hasSagePlotDir             } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { selectCurrentOrExisting    } from '../utils_nfcore_oncoanalyser_pipeline/utils'

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

            def redux_dir_tumor_selected = selectCurrentOrExisting(redux_dir_tumor, getTumorDnaReduxDir(meta))
            def redux_dir_normal_selected = selectCurrentOrExisting(redux_dir_normal, getNormalDnaReduxDir(meta))
            def redux_dir_donor_selected = selectCurrentOrExisting(redux_dir_donor, getDonorDnaReduxDir(meta))

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
                selectCurrentOrExisting(purple_dir, getPurpleDir(meta)),
            ]

        }
        .branch { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx, redux_tsvs, purple_dir ->

            def has_existing = hasSagePlotDir(meta)

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
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    visualiser_dir = ch_outputs  // channel: [ meta, sage_visualiser_dir ]
}
