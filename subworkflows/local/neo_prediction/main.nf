//
// Neo identifies and scores neoepitopes
//

nextflow.enable.types = true

include { NEO_ANNOTATE_FUSIONS  } from '../../../modules/local/neo/annotate_fusions/main'
include { NEO_FINDER  } from '../../../modules/local/neo/finder/main'
include { NEO_SCORER  } from '../../../modules/local/neo/scorer/main'

include { FileType                 } from '../utils_nfcore_oncoanalyser_pipeline/types'
include { groupByMeta              } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { joinMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { restoreMeta              } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { getInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSample        } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSampleName    } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorRnaSample        } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorRnaSampleName    } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { hasNormalDna             } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { hasTumorRna              } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { selectCurrentOrExisting  } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow NEO_PREDICTION {
    take:
    // Sample data
    ch_inputs                  // channel: [mandatory] [ meta ]
    ch_tumor_rna_aln           // channel: [mandatory] [ meta, aln, idx ]
    ch_isofox_dir              // channel: [mandatory] [ meta, isofox_dir ]
    ch_purple_dir              // channel: [mandatory] [ meta, purple_dir ]
    ch_sage_append_dir_somatic // channel: [mandatory] [ meta, sage_append_dir ]
    ch_lilac_dir               // channel: [mandatory] [ meta, lilac_dir ]
    ch_linx_annotation_dir     // channel: [mandatory] [ meta, linx_annotation_dir ]

    // Reference data
    genome_fasta               // channel: [mandatory] /path/to/genome_fasta
    genome_version             // channel: [mandatory] genome version
    genome_fai                 // channel: [mandatory] /path/to/genome_fai
    ensembl_data_resources     // channel: [mandatory] /path/to/ensembl_data_resources/
    neo_resources              // channel: [mandatory] /path/to/neo_resources/
    cohort_tpm_medians         // channel: [mandatory] /path/to/cohort_tpm_medians/

    // Params
    isofox_read_length         //  string: [mandatory] Isofox read length

    main:
    //
    // MODULE: Neo finder
    //
    // Select input sources
    // channel: [ meta, purple_dir, linx_annotation_dir ]
    ch_finder_inputs_selected = groupByMeta([
        ch_purple_dir,
        ch_linx_annotation_dir,
    ])
        .map { meta, purple_dir, linx_annotation_dir ->

            return [
                meta,
                selectCurrentOrExisting(purple_dir, getInput(getTumorDnaSample(meta), FileType.PURPLE_DIR)),
                selectCurrentOrExisting(linx_annotation_dir, getInput(getTumorDnaSample(meta), FileType.LINX_ANNO_DIR)),
            ]

        }

    // Sort inputs
    // channel: runnable: [ meta, purple_dir, linx_annotation_dir ]
    // channel: skip: [ meta ]
    ch_finder_inputs_sorted = ch_finder_inputs_selected
        .branch { meta, purple_dir, linx_annotation_dir ->

            def has_normal_dna = hasNormalDna(meta)

            def tumor_id = getTumorDnaSampleName(meta)
            def has_smlv_vcf = purple_dir ? purple_dir.resolve("${tumor_id}.purple.somatic.vcf.gz").exists() : false

            def has_runnable_inputs = has_smlv_vcf && linx_annotation_dir && has_normal_dna

            runnable: has_runnable_inputs
            skip: true
                return meta
        }

    // Create process input channel
    // channel: sample_data: [ meta_finder, purple_dir, linx_annotation_dir ]
    ch_finder_inputs = ch_finder_inputs_sorted.runnable
        .map { meta, purple_dir, linx_annotation_dir ->

            def meta_finder = [
                key: meta.case_id,
                id: meta.case_id,
                sample_id: getTumorDnaSampleName(meta),
            ]

            return [meta_finder, purple_dir, linx_annotation_dir]
        }

    // Run process
    NEO_FINDER(
        ch_finder_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        ensembl_data_resources,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, neo_finder_dir ]
    ch_finder_out = restoreMeta(channel.topic('neo_finder_dir'), ch_inputs)

    //
    // MODULE: Fusion annotation (Isofox)
    //
    // Annotate the fusion-derived neoepitope using Isofox where RNA data is available

    // Select input sources and sort
    // channel: runnable: [ meta, neo_finder_dir, tumor_rna_aln, tumor_rna_idx ]
    // channel: skip: [ meta ]
    ch_isofox_inputs_sorted = groupByMeta([
        ch_finder_out,
        ch_tumor_rna_aln,
    ])
        .map { meta, neo_finder_dir, tumor_rna_aln, tumor_rna_idx ->
            return [
                meta,
                neo_finder_dir,
                selectCurrentOrExisting(tumor_rna_aln, getInput(getTumorRnaSample(meta), FileType.ALN)),
                selectCurrentOrExisting(tumor_rna_idx, getInput(getTumorRnaSample(meta), FileType.IDX)),
            ]
        }
        .branch { meta, neo_finder_dir, tumor_rna_aln, tumor_rna_idx ->
            runnable: hasTumorRna(meta)
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_isofox, neo_finder_dir, tumor_rna_aln, tumor_rna_idx ]
    ch_isofox_inputs = ch_isofox_inputs_sorted.runnable
        .map { meta, neo_finder_dir, tumor_rna_aln, tumor_rna_idx ->

            def meta_isofox = [
                key: meta.case_id,
                id: meta.case_id,
                sample_id: getTumorDnaSampleName(meta),
            ]

            return [meta_isofox, neo_finder_dir, tumor_rna_aln, tumor_rna_idx]
        }

    // Run process
    NEO_ANNOTATE_FUSIONS(
        ch_isofox_inputs,
        isofox_read_length,
        genome_fasta,
        genome_version,
        genome_fai,
        ensembl_data_resources,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, annotated_fusions ]
    ch_annotate_fusions_out = channel.empty()
        .mix(
            restoreMeta(channel.topic('neo_annotated_fusions_tsv'), ch_inputs),
            ch_isofox_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    //
    // MODULE: Neo scorer
    //
    // Select input sources and prepare input channel
    // channel: [ meta_scorer, isofox_dir, purple_dir, sage_append_dir_somatic, lilac_dir, neo_finder_dir, annotated_fusions ]
    ch_scorer_inputs = groupByMeta([
        ch_isofox_dir,
        ch_purple_dir,
        ch_sage_append_dir_somatic,
        ch_lilac_dir,
        ch_finder_out,
        ch_annotate_fusions_out,
    ])
        .map { meta, isofox_dir, purple_dir, sage_append_dir_somatic, lilac_dir, neo_finder_dir, annotated_fusions ->

            def meta_scorer = [
                key: meta.case_id,
                id: meta.case_id,
                sample_id: getTumorDnaSampleName(meta),
                cancer_type: meta.cancer_type,
            ]

            def sage_append_vcf_somatic = null
            if (hasTumorRna(meta)) {
                meta_scorer.sample_rna_id = getTumorRnaSampleName(meta)

                def sage_append_dir_somatic_selected = selectCurrentOrExisting(sage_append_dir_somatic, getInput(getTumorDnaSample(meta), FileType.SAGE_APPEND_DIR))
                sage_append_vcf_somatic = file(sage_append_dir_somatic_selected).resolve("${meta_scorer.sample_id}.sage.append.vcf.gz")
            }

            return [
                meta_scorer,
                selectCurrentOrExisting(isofox_dir, getInput(getTumorRnaSample(meta), FileType.ISOFOX_DIR)),
                selectCurrentOrExisting(purple_dir, getInput(getTumorDnaSample(meta), FileType.PURPLE_DIR)),
                sage_append_vcf_somatic,
                selectCurrentOrExisting(lilac_dir, getInput(getTumorDnaSample(meta), FileType.LILAC_DIR)),
                neo_finder_dir,
                annotated_fusions,
            ]
        }
        .branch { meta, isofox_dir, purple_dir, sage_append_dir_somatic, lilac_dir, neo_finder_dir, annotated_fusions ->
            runnable: purple_dir && neo_finder_dir && lilac_dir
            skip: true
                return meta
        }

    // Run process
    NEO_SCORER(
        ch_scorer_inputs.runnable,
        ensembl_data_resources,
        neo_resources,
        cohort_tpm_medians,
    )

    // channel: [ meta, neo_finder_dir ]
    // channel: [ meta, annotated_fusions ]
    // channel: [ meta, neo_scorer_dir ]
    emit:
    finder_out = ch_finder_out
    annotated_fusions_out = ch_annotate_fusions_out
    neo_scorer_dir = restoreMeta(channel.topic('neo_scorer_dir'), ch_inputs)
}
