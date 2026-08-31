//
// CUPPA predicts tissue of origin from molecular profiles
//

nextflow.enable.types = true

include { CUPPA } from '../../../modules/local/cuppa/main'

include { getPurpleSomaticVcf     } from '../utils_nfcore_oncoanalyser_pipeline/accessors_outputs'
include { getInput                } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSample       } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSampleName   } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorRnaSample       } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorRnaSampleName   } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasInput                } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasNormalDna            } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasTumorDna             } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasTumorRna             } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { FileType                } from '../utils_nfcore_oncoanalyser_pipeline/types_enums'
include { groupByMeta             } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { joinMeta                } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { restoreMeta             } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { selectCurrentOrExisting } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow CUPPA_PREDICTION {
    take:
    // Sample data
    ch_inputs              : Channel<Map>              // channel: [mandatory] [ meta ]
    ch_isofox_dir          : Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, isofox_dir ]
    ch_purple_dir          : Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, purple_dir ]
    ch_linx_annotation_dir : Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, linx_annotation_dir ]
    ch_virusinterpreter_dir: Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, virusinterpreter_dir ]

    // Reference data
    genome_version         : Channel<String>           // channel: [mandatory] genome version
    cuppa_alt_sj           : Channel<Path>             // channel: [mandatory] /path/to/cuppa_alt_sj/
    cuppa_classifier       : Channel<Path>             // channel: [mandatory] /path/to/cuppa_classifier/

    main:
    // Select input sources then sort
    // channel: runnable: [ meta, isofox_dir, purple_dir, linx_annotation_dir, virusinterpreter_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = groupByMeta([
        ch_isofox_dir,
        ch_purple_dir,
        ch_linx_annotation_dir,
        ch_virusinterpreter_dir,
    ])
        .map { meta, isofox_dir, purple_dir, linx_annotation_dir, virusinterpreter_dir ->
            return [
                meta,
                selectCurrentOrExisting(isofox_dir, getInput(getTumorRnaSample(meta), FileType.ISOFOX_DIR)),
                selectCurrentOrExisting(purple_dir, getInput(getTumorDnaSample(meta), FileType.PURPLE_DIR)),
                selectCurrentOrExisting(linx_annotation_dir, getInput(getTumorDnaSample(meta), FileType.LINX_ANNO_DIR)),
                selectCurrentOrExisting(virusinterpreter_dir, getInput(getTumorDnaSample(meta), FileType.VIRUSINTERPRETER_DIR)),
            ]
        }
        .branch { meta, isofox_dir, purple_dir, linx_annotation_dir, virusinterpreter_dir ->

            // Run the following:
            //   - tumor DNA and normal DNA
            //   - tumor DNA and normal DNA, and tumor RNA
            //   - tumor RNA only
            //
            // Do not run the following:
            //   - tumor DNA only
            //   - panel mode (controlled by excluded from targeted subworkflow)
            //
            // (run exclusions currently done basis for presence of normal DNA)

            def has_existing = hasInput(getTumorDnaSample(meta), FileType.CUPPA_DIR)
            def has_normal_dna = hasNormalDna(meta)

            def tumor_dna_id = getTumorDnaSampleName(meta)
            def has_smlv_vcf = getPurpleSomaticVcf(tumor_dna_id, purple_dir)?.exists() ?: false

            def has_runnable_inputs = isofox_dir || (has_smlv_vcf && linx_annotation_dir && has_normal_dna)

            runnable: has_runnable_inputs && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: sample_data: [ meta_cuppa, isofox_dir, purple_dir, linx_annotation_dir, virusinterpreter_dir ]
    // channel: categories: [ categories ]
    ch_cuppa_inputs = ch_inputs_sorted.runnable
        .multiMap { meta, isofox_dir, purple_dir, linx_annotation_dir, virusinterpreter_dir ->

            def has_tumor_dna = hasTumorDna(meta)
            def has_normal_dna = hasNormalDna(meta)
            def has_tumor_rna = hasTumorRna(meta)

            def has_dna_inputs = (purple_dir && linx_annotation_dir)
            def has_rna_inputs = isofox_dir

            def run_dna = has_dna_inputs && has_tumor_dna && has_normal_dna
            def run_rna = has_rna_inputs && has_tumor_rna

            def tumor_dna_id = getTumorDnaSampleName(meta)
            def tumor_rna_id = getTumorRnaSampleName(meta)

            def categories
            def sample_id = null
            def sample_rna_id = null

            if (run_dna && run_rna) {

                categories = 'ALL'

                sample_id = tumor_dna_id
                sample_rna_id = tumor_rna_id

            } else if (run_dna) {

                categories = 'DNA'

                sample_id = tumor_dna_id

                isofox_dir = null

            } else if (run_rna) {

                categories = 'RNA'

                sample_id = has_tumor_dna ? tumor_dna_id : tumor_rna_id
                sample_rna_id = tumor_rna_id

                purple_dir = null
                linx_annotation_dir = null
                virusinterpreter_dir = null

            } else {

                assert false

            }

            def meta_cuppa = record(
                key: meta.case_id,
                id: meta.case_id,
                sample_id: sample_id,
                sample_rna_id: sample_rna_id,
            )

            sample_data: [meta_cuppa, isofox_dir, purple_dir, linx_annotation_dir, virusinterpreter_dir]
            categories: categories
    }

    // Run process
    CUPPA(
        ch_cuppa_inputs.sample_data,
        genome_version,
        cuppa_alt_sj,
        cuppa_classifier,
        ch_cuppa_inputs.categories,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, cuppa_dir ]
    ch_outputs = channel.empty()
        .mix(
            restoreMeta(channel.topic('cuppa_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    emit:
    cuppa_dir = ch_outputs // channel: [ meta, cuppa_dir ]
}
