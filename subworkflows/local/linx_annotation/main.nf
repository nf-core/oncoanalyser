//
// LINX annotates and interprets structural variants
//

include { LINX_GERMLINE } from '../../../modules/local/linx/germline/main'
include { LINX_SOMATIC  } from '../../../modules/local/linx/somatic/main'

workflow LINX_ANNOTATION {
    take:
    // Sample data
    ch_inputs              // channel: [mandatory] [ meta ]
    ch_purple_dir          // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_version         // channel: [mandatory] genome version
    ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/
    known_fusion_data      // channel: [mandatory] /path/to/known_fusion_data
    driver_gene_panel      // channel: [mandatory] /path/to/driver_gene_panel

    main:
    //
    // STEP: Handle inputs
    //
    // Select input sources then sort
    // channel: runnable: [ meta, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_purple_dir
        .map { meta, purple_dir ->
            return [
                meta,
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
            ]
        }
        .branch { meta, purple_dir ->
            runnable: purple_dir
            skip: true
                return meta
        }

    //
    // MODULE: LINX germline annotation
    //
    // Select inputs that are eligible to run
    // channel: runnable: [ meta, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_germline_sorted = ch_inputs_sorted.runnable
        .branch { meta, purple_dir ->

            def tumor_id = Utils.getTumorDnaSampleName(meta)

            def has_tumor_normal = Utils.hasTumorDna(meta) && Utils.hasNormalDna(meta)
            def has_sv_germline_vcf = purple_dir.resolve("${tumor_id}.purple.sv.germline.vcf.gz").exists()
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.LINX_ANNO_DIR_NORMAL)

            runnable: has_tumor_normal && has_sv_germline_vcf && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta, purple_sv_vcf ]
    ch_linx_germline_inputs = ch_inputs_germline_sorted.runnable
        .map { meta, purple_dir ->

            def tumor_id = Utils.getTumorDnaSampleName(meta)

            def meta_linx = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: tumor_id,
            ]

            def sv_vcf = purple_dir.resolve("${tumor_id}.purple.sv.germline.vcf.gz")

            return [meta_linx, sv_vcf]
        }

    // Run process
    LINX_GERMLINE(
        ch_linx_germline_inputs,
        genome_version,
        ensembl_data_resources,
        driver_gene_panel,
    )

    //
    // MODULE: LINX somatic annotation
    //
    // Select inputs that are eligible to run
    // channel: runnable: [ meta, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_somatic_sorted = ch_inputs_sorted.runnable
        .branch { meta, purple_dir ->

            def has_tumor = Utils.hasTumorDna(meta)
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.LINX_ANNO_DIR_TUMOR)

            runnable: has_tumor && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta, purple_dir ]
    ch_linx_somatic_inputs = ch_inputs_somatic_sorted.runnable
        .map { meta, purple_dir ->

            def meta_linx = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
            ]

            return [meta_linx, purple_dir]
        }

    // Run process
    LINX_SOMATIC(
        ch_linx_somatic_inputs,
        genome_version,
        ensembl_data_resources,
        known_fusion_data,
        driver_gene_panel,
    )

    //
    // STEP: Handle outputs
    //
    // Set outputs, restoring original meta
    // channel: [ meta, linx_annotation_dir ]
    ch_outputs_somatic = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('linx_somatic_annotation_dir'), ch_inputs),
            ch_inputs_somatic_sorted.skip.map { meta -> [meta, []] },
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    ch_outputs_germline = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('linx_germline_annotation_dir'), ch_inputs),
            ch_inputs_germline_sorted.skip.map { meta -> [meta, []] },
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    somatic_dir  = ch_outputs_somatic  // channel: [ meta, linx_annotation_dir ]
    germline_dir = ch_outputs_germline // channel: [ meta, linx_annotation_dir ]
}
