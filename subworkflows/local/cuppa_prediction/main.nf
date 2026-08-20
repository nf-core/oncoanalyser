//
// CUPPA predicts tissue of origin from molecular profiles
//

include { CUPPA } from '../../../modules/local/cuppa/main'
include { groupByMeta; joinMeta; restoreMeta } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { getPurpleDir; getTumorDnaLinxAnnoDir; getTumorDnaSampleName; getTumorRnaIsofoxDir; getTumorRnaSampleName; getVirusinterpreterDir; hasCuppaDir; hasNormalDna; hasTumorDna; hasTumorRna; selectCurrentOrExisting } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow CUPPA_PREDICTION {
    take:
    // Sample data
    ch_inputs               // channel: [mandatory] [ meta ]
    ch_isofox_dir           // channel: [mandatory] [ meta, isofox_dir ]
    ch_purple_dir           // channel: [mandatory] [ meta, purple_dir ]
    ch_linx_annotation_dir  // channel: [mandatory] [ meta, linx_annotation_dir ]
    ch_virusinterpreter_dir // channel: [mandatory] [ meta, virusinterpreter_dir ]

    // Reference data
    genome_version          // channel: [mandatory] genome version
    cuppa_alt_sj            // channel: [mandatory] /path/to/cuppa_alt_sj/
    cuppa_classifier        // channel: [mandatory] /path/to/cuppa_classifier/

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
                selectCurrentOrExisting(isofox_dir, getTumorRnaIsofoxDir(meta)),
                selectCurrentOrExisting(purple_dir, getPurpleDir(meta)),
                selectCurrentOrExisting(linx_annotation_dir, getTumorDnaLinxAnnoDir(meta)),
                selectCurrentOrExisting(virusinterpreter_dir, getVirusinterpreterDir(meta)),
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

            def has_existing = hasCuppaDir(meta)
            def has_normal_dna = hasNormalDna(meta)

            def tumor_dna_id = getTumorDnaSampleName(meta)
            def has_smlv_vcf = purple_dir ? purple_dir.resolve("${tumor_dna_id}.purple.somatic.vcf.gz").exists() : false

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

            def meta_cuppa = [
                key: meta.case_id,
                id: meta.case_id,
            ]

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

            if (run_dna && run_rna) {

                categories = 'ALL'

                meta_cuppa.sample_id = tumor_dna_id
                meta_cuppa.sample_rna_id = tumor_rna_id

            } else if (run_dna) {

                categories = 'DNA'

                meta_cuppa.sample_id = tumor_dna_id

                isofox_dir = []

            } else if (run_rna) {

                categories = 'RNA'

                meta_cuppa.sample_id = has_tumor_dna ? tumor_dna_id : tumor_rna_id
                meta_cuppa.sample_rna_id = tumor_rna_id

                purple_dir = []
                linx_annotation_dir = []
                virusinterpreter_dir = []

            } else {

                assert false

            }

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
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    cuppa_dir = ch_outputs // channel: [ meta, cuppa_dir ]
}
