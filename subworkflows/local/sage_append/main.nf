//
// SAGE append adds additional sample data to an existing SAGE VCF
//

nextflow.enable.types = true

include { SAGE_APPEND as SAGE_APPEND_SOMATIC  } from '../../../modules/local/sage/append/main'
include { SAGE_APPEND as SAGE_APPEND_GERMLINE } from '../../../modules/local/sage/append/main'


include { getTumorReduxDirAlignment } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { getTumorReduxTsvs         } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { getPurpleGermlineVcf      } from '../utils_nfcore_oncoanalyser_pipeline/accessors_outputs'
include { getPurpleSomaticVcf       } from '../utils_nfcore_oncoanalyser_pipeline/accessors_outputs'
include { getInput                  } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getLongitudinalSampleName } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalDnaSample        } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalDnaSampleName    } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSample         } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSampleName     } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorRnaSample         } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorRnaSampleName     } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasInput                  } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasNormalDna              } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasTumorDna               } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasTumorRna               } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { groupByMeta               } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { joinMeta                  } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { restoreMeta               } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { FileType                  } from '../utils_nfcore_oncoanalyser_pipeline/types_enums'
include { selectCurrentOrExisting   } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow SAGE_APPEND {
    take:
    // Sample data
    ch_inputs: Channel<Map>           // channel: [mandatory] [ meta ]
    ch_purple_dir: Channel<Tuple<Map, Path>>       // channel: [mandatory] [ meta, purple_dir ]
    ch_redux_dir_tumor: Channel<Tuple<Map, Path>>  // channel: [mandatory] [ meta, redux_dir ]
    ch_tumor_rna_aln: Channel<Tuple<Map, Path, Path>>    // channel: [mandatory] [ meta, aln, idx ]

    // Reference data
    genome_fasta: Channel<Path>        // channel: [mandatory] /path/to/genome_fasta
    genome_version: Channel<String>      // channel: [mandatory] genome version
    genome_fai: Channel<Path>          // channel: [mandatory] /path/to/genome_fai
    genome_dict: Channel<Path>         // channel: [mandatory] /path/to/genome_dict

    // Params
    sequencing_platform: String // string:  [mandatory] sequencing platform
    enable_germline: Boolean     // boolean: [mandatory] Enable germline
    targeted_mode: Boolean       // boolean: [mandatory] Set targeted mode
    purity_estimate_mode: Boolean // boolean: [mandatory] Set purity estimate mode

    main:
    //
    // STEP: Handle inputs
    //

    // Select input sources then sort
    // channel: runnable: [ meta, purple_dir, tumor_dna_aln, tumor_dna_idx, [redux_tsv_tumor, ...], tumor_rna_aln, tumor_rna_idx ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = groupByMeta([
        ch_purple_dir,
        ch_redux_dir_tumor,
        ch_tumor_rna_aln,
    ])
        .map { meta, purple_dir, redux_dir_tumor, tumor_rna_aln, tumor_rna_idx ->

            def redux_dir_tumor_selected = selectCurrentOrExisting(redux_dir_tumor, getInput(getTumorDnaSample(meta), FileType.REDUX_DIR))
            def (tumor_aln, tumor_idx) = getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)
            def redux_tsvs_tumor = getTumorReduxTsvs(meta, redux_dir_tumor_selected)

            return [
                meta,
                selectCurrentOrExisting(purple_dir, getInput(getTumorDnaSample(meta), FileType.PURPLE_DIR)),
                tumor_aln,
                tumor_idx,
                redux_tsvs_tumor,
                selectCurrentOrExisting(tumor_rna_aln, getInput(getTumorRnaSample(meta), FileType.ALN)),
                selectCurrentOrExisting(tumor_rna_idx, getInput(getTumorRnaSample(meta), FileType.IDX)),
            ]

        }
        .branch { meta, purple_dir, tumor_dna_aln, tumor_dna_idx, redux_tsvs_tumor, tumor_rna_aln, tumor_rna_idx ->
            def has_aln = tumor_dna_aln || tumor_rna_aln
            runnable: has_aln && purple_dir
            skip: true
                return meta
        }

    //
    // MODULE: SAGE append germline
    //
    // Select inputs that are eligible to run
    // channel: runnable: [ meta, purple_dir, tumor_dna_aln, tumor_dna_idx, [redux_tsv_tumor, ...], tumor_rna_aln, tumor_rna_idx ]
    // channel: skip: [ meta ]
    ch_inputs_germline_sorted = ch_inputs_sorted.runnable
        .branch { meta, purple_dir, tumor_dna_aln, tumor_dna_idx, redux_tsvs_tumor, tumor_rna_aln, tumor_rna_idx ->

            // NOTE(SW): explicit in expectation to always obtain the primary tumor DNA sample ID here
            def tumor_dna_id = getTumorDnaSampleName(meta)

            def has_tumor_rna = hasTumorRna(meta)
            def has_normal_dna = hasNormalDna(meta)
            def has_smlv_germline = getPurpleGermlineVcf(tumor_dna_id, purple_dir).exists()

            def should_append_rna_variants = has_tumor_rna && has_normal_dna && has_smlv_germline

            def has_existing = hasInput(getNormalDnaSample(meta), FileType.SAGE_APPEND_DIR)

            runnable: should_append_rna_variants && ! has_existing && enable_germline
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_append, purple_smlv_vcf, [aln, ...], [idx, ...], [redux_tsvs_tumor, ...] ]
    ch_sage_append_germline_inputs = ch_inputs_germline_sorted.runnable
        .map { meta, purple_dir, _tumor_dna_aln, _tumor_dna_idx, _redux_tsvs_tumor, tumor_rna_aln, tumor_rna_idx ->

            // NOTE(SW): explicit in expectation to always obtain the primary tumor DNA sample ID here
            def tumor_dna_id = getTumorDnaSampleName(meta)
            def output_file_id = getNormalDnaSampleName(meta)

            def meta_append = [
                key: meta.case_id,
                topic_key: 'germline',
                id: "${meta.case_id}:${output_file_id}",
                output_file_id: output_file_id,
                reference_ids: [getTumorRnaSampleName(meta)],
            ]

            def alns = [tumor_rna_aln]
            def idxs = [tumor_rna_idx]

            def purple_smlv_vcf = getPurpleGermlineVcf(tumor_dna_id, purple_dir)

            // NOTE(SW): in stub mode we allow absence of indexes so must handle here
            idxs = idxs.flatten()

            return [meta_append, purple_smlv_vcf, alns, idxs, []]

        }

    // Run process
    SAGE_APPEND_GERMLINE(
        ch_sage_append_germline_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        sequencing_platform,
        targeted_mode,
    )

    //
    // MODULE: SAGE append somatic
    //
    // Select inputs that are eligible to run
    // channel: runnable: [ meta, purple_dir, tumor_dna_aln, tumor_dna_idx, [redux_tsv_tumor, ...], tumor_rna_aln, tumor_rna_idx ]
    // channel: skip: [ meta ]
    ch_inputs_somatic_sorted = ch_inputs_sorted.runnable
        .branch { meta, purple_dir, tumor_dna_aln, tumor_dna_idx, redux_tsvs_tumor, tumor_rna_aln, tumor_rna_idx ->

            def tumor_dna_id = getTumorDnaSampleName(meta)

            def has_tumor_rna = hasTumorRna(meta)
            def has_tumor_dna = hasTumorDna(meta)
            def has_smlv_somatic = getPurpleSomaticVcf(tumor_dna_id, purple_dir).exists()

            def should_append_rna_variants = ! purity_estimate_mode && has_tumor_rna && has_tumor_dna && has_smlv_somatic
            def should_append_longitudinal_variants = purity_estimate_mode && has_tumor_dna && has_smlv_somatic

            def has_existing = hasInput(getTumorDnaSample(meta), FileType.SAGE_APPEND_DIR)

            runnable: (should_append_rna_variants || should_append_longitudinal_variants) && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_append, purple_smlv_vcf, [aln, ...], [idx, ...], [redux_tsvs_tumor, ...] ]
    ch_sage_append_somatic_inputs = ch_inputs_somatic_sorted.runnable
        .map { meta, purple_dir, tumor_dna_aln, tumor_dna_idx, redux_tsvs_tumor, tumor_rna_aln, tumor_rna_idx ->

            // NOTE(SW): explicit in expectation to always obtain the primary tumor DNA sample ID here
            def tumor_dna_id = getTumorDnaSampleName(meta)
            def output_file_id = purity_estimate_mode ? getLongitudinalSampleName(meta) : tumor_dna_id

            def meta_append = [
                key: meta.case_id,
                topic_key: 'somatic',
                id: "${meta.case_id}:${output_file_id}",
                output_file_id: output_file_id,
                reference_ids: [],
            ]

            def alns = []
            def idxs = []
            def redux_tsvs = []

            if (! purity_estimate_mode && tumor_rna_aln) {
                meta_append.reference_ids.add(getTumorRnaSampleName(meta))
                alns.add(tumor_rna_aln)
                idxs.add(tumor_rna_idx)
            }

            if (purity_estimate_mode && tumor_dna_aln) {
                meta_append.reference_ids.add(getLongitudinalSampleName(meta))
                alns.add(tumor_dna_aln)
                idxs.add(tumor_dna_idx)
                redux_tsvs = redux_tsvs_tumor
            }

            def purple_smlv_vcf = getPurpleSomaticVcf(tumor_dna_id, purple_dir)

            // NOTE(SW): in stub mode we allow absence of indexes so must handle here
            idxs = idxs.flatten()

            return [meta_append, purple_smlv_vcf, alns, idxs, redux_tsvs]
        }

    // Run process
    SAGE_APPEND_SOMATIC(
        ch_sage_append_somatic_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        sequencing_platform,
        targeted_mode,
    )

    //
    // STEP: Handle outputs
    //
    // Set outputs, restoring original meta
    // NOTE(SW): subscribe once and reuse, channel.topic() is single-consumer under types
    ch_sage_append_topic = channel.topic('sage_append_dir')

    // channel: [ meta, sage_append_dir ]
    ch_outputs_somatic = channel.empty()
        .mix(
            restoreMeta(ch_sage_append_topic.filter { d -> d[0].topic_key == 'somatic' }, ch_inputs),
            ch_inputs_somatic_sorted.skip.map { meta -> [meta, null] },
            ch_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    // channel: [ meta, sage_append_dir ]
    ch_outputs_germline = channel.empty()
        .mix(
            restoreMeta(ch_sage_append_topic.filter { d -> d[0].topic_key == 'germline' }, ch_inputs),
            ch_inputs_germline_sorted.skip.map { meta -> [meta, null] },
            ch_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    emit:
    somatic_dir  = ch_outputs_somatic  // channel: [ meta, sage_append_dir ]
    germline_dir = ch_outputs_germline // channel: [ meta, sage_append_dir ]
}
