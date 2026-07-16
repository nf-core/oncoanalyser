//
// CUPPA predicts tissue of origin from molecular profiles
//

include { CUPPA } from '../../../modules/local/cuppa/main'

workflow CUPPA_PREDICTION {
    take:
    // Sample data
    ch_inputs           // channel: [mandatory] [ meta ]
    ch_isofox           // channel: [mandatory] [ meta, isofox_dir ]
    ch_purple           // channel: [mandatory] [ meta, purple_dir ]
    ch_linx             // channel: [mandatory] [ meta, linx_annotation_dir ]
    ch_virusinterpreter // channel: [mandatory] [ meta, virusinterpreter_dir ]

    // Reference data
    genome_version      // channel: [mandatory] genome version
    cuppa_alt_sj        // channel: [mandatory] /path/to/cuppa_alt_sj/
    cuppa_classifier    // channel: [mandatory] /path/to/cuppa_classifier/

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources
    // channel: [ meta, isofox_dir, purple_dir, linx_annotation_dir, virusinterpreter_dir ]
    ch_inputs_selected = channels.WorkflowChannels.groupByMeta(
        ch_isofox,
        ch_purple,
        ch_linx,
        ch_virusinterpreter,
    )
        .map { meta, isofox_dir, purple_dir, linx_annotation_dir, virusinterpreter_dir ->

            def inputs = [
                sample.Inputs.preferUserProvidedInput(isofox_dir, meta, sample.FileKey.ISOFOX_DIR),
                sample.Inputs.preferUserProvidedInput(purple_dir, meta, sample.FileKey.PURPLE_DIR),
                sample.Inputs.preferUserProvidedInput(linx_annotation_dir, meta, sample.FileKey.LINX_ANNO_DIR_TUMOR),
                sample.Inputs.preferUserProvidedInput(virusinterpreter_dir, meta, sample.FileKey.VIRUSINTERPRETER_DIR),
            ]

            return [meta, *inputs]
        }

    // Sort inputs
    // channel: runnable: [ meta, isofox_dir, purple_dir, linx_annotation_dir, virusinterpreter_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_inputs_selected
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

            def has_existing = sample.Inputs.hasExisting(meta, sample.FileKey.CUPPA_DIR)
            def has_normal_dna = sample.Inputs.hasNormalDna(meta)

            def has_runnable_inputs = isofox_dir || (purple_dir && linx_annotation_dir && has_normal_dna)

            runnable: has_runnable_inputs && !has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: sample_data: [ meta, isofox_dir, purple_dir, linx_annotation_dir, virusinterpreter_dir ]
    // channel: categories: [ categories ]
    ch_cuppa_inputs = ch_inputs_sorted.runnable
        .multiMap{ meta, isofox_dir, purple_dir, linx_annotation_dir, virusinterpreter_dir ->

            def meta_cuppa = [
                key: meta.group_id,
                id: meta.group_id,
            ]

            def tumor_dna_id = sample.Inputs.getTumorDnaSampleName(meta)
            def tumor_rna_id = sample.Inputs.getTumorRnaSampleOutputId(meta)

            meta_cuppa.sample_id = tumor_dna_id ?: tumor_rna_id
            meta_cuppa.sample_rna_id = tumor_rna_id

            def has_tumor_dna = sample.Inputs.hasTumorDna(meta)
            def has_normal_dna = sample.Inputs.hasNormalDna(meta)
            def has_tumor_rna = sample.Inputs.hasTumorRna(meta)

            def has_dna_inputs = (purple_dir && linx_annotation_dir)
            def has_rna_inputs = isofox_dir

            def run_dna = has_dna_inputs && has_tumor_dna && has_normal_dna
            def run_rna = has_rna_inputs && has_tumor_rna

            def categories

            if (run_dna && run_rna) {
                categories = 'ALL'
            } else if (run_dna) {
                categories = 'DNA'
            } else if (run_rna) {
                categories = 'RNA'
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

    ch_versions = ch_versions.mix(CUPPA.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, cuppa_dir ]
    ch_outputs = Channel.empty()
        .mix(
            channels.WorkflowChannels.restoreMeta(CUPPA.out.cuppa_dir, ch_inputs),
            channels.PlaceholderChannels.toolDir(ch_inputs_sorted.skip),
        )

    emit:
    cuppa_dir = ch_outputs  // channel: [ meta, cuppa_dir ]

    versions  = ch_versions // channel: [ versions.yml ]
}
