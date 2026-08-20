//
// LILAC is a WGS tool for HLA typing and somatic CNV and SNV calling
//

include { LILAC } from '../../../modules/local/lilac/main'

include { FileType } from '../utils_nfcore_oncoanalyser_pipeline/types'
include { groupByMeta                } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { joinMeta                   } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { restoreMeta                } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { getInput                   } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getNormalDnaSample         } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getNormalDnaSampleName     } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getNormalReduxDirAlignment } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getTumorDnaSample          } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getTumorDnaSampleName      } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getTumorReduxDirAlignment  } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getTumorRnaSample          } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { hasInput                   } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { hasNormalDna               } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { hasTumorDna                } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { selectCurrentOrExisting    } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow LILAC_CALLING {
    take:
    // Sample data
    ch_inputs           // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor  // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal // channel: [mandatory] [ meta, redux_dir ]
    ch_tumor_rna_aln    // channel: [mandatory] [ meta, aln, idx ]
    ch_purple           // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_fasta        // channel: [mandatory] /path/to/genome_fasta
    genome_version      // channel: [mandatory] genome version
    genome_fai          // channel: [mandatory] /path/to/genome_fai
    lilac_resource_dir  // channel: [mandatory] /path/to/lilac_resource_dir/

    // Params
    sequencing_platform // string:  [mandatory] sequencing platform
    targeted_mode       // boolean: [mandatory] Set targeted mode

    main:
    // Select input sources then sort
    // channel: runnable: [meta, normal_dna_aln, normal_dna_idx, tumor_dna_aln, tumor_dna_idx, tumor_rna_aln, tumor_rna_idx, purple_dir]
    // channel: skip: [ meta ]
    ch_dna_inputs_sorted = groupByMeta([
        ch_redux_dir_tumor,
        ch_redux_dir_normal,
        ch_tumor_rna_aln,
        ch_purple,
    ])
        .map { meta, redux_dir_tumor, redux_dir_normal, tumor_rna_aln, tumor_rna_idx, purple_dir ->

            def redux_dir_tumor_selected = selectCurrentOrExisting(redux_dir_tumor, getInput(getTumorDnaSample(meta), FileType.REDUX_DIR))
            def redux_dir_normal_selected = selectCurrentOrExisting(redux_dir_normal, getInput(getNormalDnaSample(meta), FileType.REDUX_DIR))

            def (tumor_dna_aln, tumor_dna_idx) = getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)
            def (normal_dna_aln, normal_dna_idx) = getNormalReduxDirAlignment(meta, redux_dir_normal_selected)

            return [
                meta,
                normal_dna_aln,
                normal_dna_idx,
                tumor_dna_aln,
                tumor_dna_idx,
                selectCurrentOrExisting(tumor_rna_aln, getInput(getTumorRnaSample(meta), FileType.ALN)),
                selectCurrentOrExisting(tumor_rna_idx, getInput(getTumorRnaSample(meta), FileType.IDX)),
                selectCurrentOrExisting(purple_dir, getInput(getTumorDnaSample(meta), FileType.PURPLE_DIR)),
            ]

        }
        .branch { meta, normal_dna_aln, normal_dna_idx, tumor_dna_aln, tumor_dna_idx, tumor_rna_aln, tumor_rna_idx, purple_dir ->

            def has_existing = hasInput(getTumorDnaSample(meta), FileType.LILAC_DIR)

            def tumor_normal_mode = tumor_dna_aln && normal_dna_aln

            def tumor_dna_id = getTumorDnaSampleName(meta)
            def has_tn_smlv_vcf = purple_dir ? purple_dir.resolve("${tumor_dna_id}.purple.somatic.vcf.gz").exists() : false

            runnable: (tumor_dna_aln || normal_dna_aln) && (has_tn_smlv_vcf || ! tumor_normal_mode) && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_lilac, normal_dna_aln, normal_dna_idx, tumor_dna_aln, tumor_dna_idx, tumor_rna_aln, tumor_rna_idx, purple_dir
    ch_lilac_inputs = groupByMeta([
        ch_dna_inputs_sorted.runnable,
    ])
        .map { meta, normal_dna_aln, normal_dna_idx, tumor_dna_aln, tumor_dna_idx, tumor_rna_aln, tumor_rna_idx, purple_dir ->

            def meta_lilac = [
                key: meta.case_id,
                id: meta.case_id,
            ]

            if (hasTumorDna(meta)) {
                meta_lilac.tumor_id = getTumorDnaSampleName(meta)
            }

            if (hasNormalDna(meta)) {
                meta_lilac.normal_id = getNormalDnaSampleName(meta)
            }

            return [meta_lilac, normal_dna_aln, normal_dna_idx, tumor_dna_aln, tumor_dna_idx, tumor_rna_aln, tumor_rna_idx, purple_dir]
        }

    // Run process
    LILAC(
        ch_lilac_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        lilac_resource_dir,
        targeted_mode,
        sequencing_platform,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, lilac_dir ]
    ch_outputs = channel.empty()
        .mix(
            restoreMeta(channel.topic('lilac_dir'), ch_inputs),
            ch_dna_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    lilac_dir = ch_outputs // channel: [ meta, lilac_dir ]
}
