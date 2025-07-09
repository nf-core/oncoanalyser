//
// LINX annotates and interprets structural variants
//

import Constants
import Utils

include { LINX_GERMLINE } from '../../../modules/local/linx/germline/main'
include { LINX_SOMATIC  } from '../../../modules/local/linx/somatic/main'

workflow LINX_ANNOTATION {
    take:
    // Sample data
    ch_inputs              // channel: [mandatory] [ meta ]
    ch_purple              // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_version         // channel: [mandatory] genome version
    ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/
    known_fusion_data      // channel: [mandatory] /path/to/known_fusion_data
    driver_gene_panel      // channel: [mandatory] /path/to/driver_gene_panel

    main:
    // Channel for versions.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources and sort
    // channel: runnable: [ meta, purple_dir ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_purple
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
            def has_sv_germline_vcf = file(purple_dir).resolve("${tumor_id}.purple.sv.germline.vcf.gz")
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.LINX_ANNO_DIR_NORMAL)

            runnable: has_tumor_normal && has_sv_germline_vcf && !has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta, sv_vcf ]
    ch_linx_germline_inputs = ch_inputs_germline_sorted.runnable
        .map { meta, purple_dir ->

            def tumor_id = Utils.getTumorDnaSampleName(meta)

            def meta_linx = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: tumor_id,
            ]

            def sv_vcf = file(purple_dir).resolve("${tumor_id}.purple.sv.germline.vcf.gz")

            return [meta_linx, sv_vcf]
        }

    // Run process
    LINX_GERMLINE(
        ch_linx_germline_inputs,
        genome_version,
        ensembl_data_resources,
        driver_gene_panel,
    )

    ch_versions = ch_versions.mix(LINX_GERMLINE.out.versions)

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

            runnable: has_tumor && !has_existing
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

    ch_versions = ch_versions.mix(LINX_SOMATIC.out.versions)


    // Set outputs, restoring original meta
    // channel: [ meta, linx_annotation_dir ]
    ch_somatic_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(LINX_SOMATIC.out.annotation_dir, ch_inputs),
            ch_inputs_somatic_sorted.skip.map { meta -> [meta, []] },
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    ch_germline_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(LINX_GERMLINE.out.annotation_dir, ch_inputs),
            ch_inputs_germline_sorted.skip.map { meta -> [meta, []] },
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    somatic  = ch_somatic_out  // channel: [ meta, linx_annotation_dir ]
    germline = ch_germline_out // channel: [ meta, linx_annotation_dir ]

    versions = ch_versions     // channel: [ versions.yml ]
}
