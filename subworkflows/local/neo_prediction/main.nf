//
// Neo identifies and scores neoepitopes
//

import Constants
import Utils

include { NEO_ANNOTATE_FUSIONS } from '../../../modules/local/neo/annotate_fusions/main'
include { NEO_FINDER           } from '../../../modules/local/neo/finder/main'
include { NEO_SCORER           } from '../../../modules/local/neo/scorer/main'

workflow NEO_PREDICTION {
    take:
    // Sample data
    ch_inputs              // channel: [mandatory] [ meta ]
    ch_tumor_rna_bam       // channel: [mandatory] [ meta, bam, bai ]
    ch_isofox              // channel: [mandatory] [ meta, isofox_dir ]
    ch_purple              // channel: [mandatory] [ meta, purple_dir ]
    ch_sage_somatic_append // channel: [mandatory] [ meta, sage_append_vcf ]
    ch_lilac               // channel: [mandatory] [ meta, lilac_dir ]
    ch_linx                // channel: [mandatory] [ meta, linx_annotation_dir ]

    // Reference data
    genome_fasta           // channel: [mandatory] /path/to/genome_fasta
    genome_version         // channel: [mandatory] genome version
    genome_fai             // channel: [mandatory] /path/to/genome_fai
    ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/
    neo_resources          // channel: [mandatory] /path/to/neo_resources/
    cohort_tpm_medians     // channel: [mandatory] /path/to/cohort_tpm_medians/

    // Params
    isofox_read_length     //  string: [mandatory] Isofox read length

    main:
    // Channel for versions.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    //
    // MODULE: Neo finder
    //
    // Select input sources
    // channel: [ meta, purple_dir, linx_annotation_dir ]
    ch_finder_inputs_selected = WorkflowOncoanalyser.groupByMeta(
        ch_purple,
        ch_linx,
    )
        .map { meta, purple_dir, linx_annotation_dir ->

            def inputs = [
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
                Utils.selectCurrentOrExisting(linx_annotation_dir, meta, Constants.INPUT.LINX_ANNO_DIR_TUMOR),
            ]

            return [meta, *inputs]
        }

    // Sort inputs
    // channel: runnable: [ meta, purple_dir, linx_annotation_dir ]
    // channel: skip: [ meta ]
    ch_finder_inputs_sorted = ch_finder_inputs_selected
        .branch { meta, purple_dir, linx_annotation_dir ->

            def has_normal_dna = Utils.hasNormalDna(meta)

            def has_runnable_inputs = purple_dir && linx_annotation_dir && has_normal_dna

            runnable: has_runnable_inputs
            skip: true
                return meta
        }

    // Create process input channel
    // channel: sample_data: [ meta_finder, purple_dir, linx_annotation_dir ]
    ch_finder_inputs = ch_finder_inputs_sorted.runnable
        .map { meta, purple_dir, linx_annotation_dir ->

            def meta_finder = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
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

    ch_versions = ch_versions.mix(NEO_FINDER.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, neo_finder_dir ]
    ch_finder_out = WorkflowOncoanalyser.restoreMeta(NEO_FINDER.out.neo_finder_dir, ch_inputs)

    //
    // MODULE: Fusion annotation (Isofox)
    //
    // Annotate the fusion-derived neoepitope using Isofox where RNA data is available

    // Select input sources and sort
    // channel: runnable: [ meta, neo_finder_dir, tumor_bam_rna, tumor_bai_rna ]
    // channel: skip: [ meta ]
    ch_isofox_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_finder_out,
        ch_tumor_rna_bam,
    )
        .map { meta, neo_finder_dir, tumor_bam, tumor_bai ->
            return [
                meta,
                neo_finder_dir,
                Utils.selectCurrentOrExisting(tumor_bam, meta, Constants.INPUT.BAM_RNA_TUMOR),
                Utils.selectCurrentOrExisting(tumor_bai, meta, Constants.INPUT.BAI_RNA_TUMOR),
            ]
        }
        .branch { meta, neo_finder_dir, tumor_bam, tumor_bai ->
            runnable: Utils.hasTumorRna(meta)
                return [meta, neo_finder_dir, tumor_bam, tumor_bai]
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_isofox, neo_finder_dir, tumor_bam_rna, tumor_bai_rna ]
    ch_isofox_inputs = ch_isofox_inputs_sorted.runnable
        .map { meta, neo_finder_dir, tumor_bam_rna, tumor_bai_rna ->

            def meta_isofox = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
            ]

            return [meta_isofox, neo_finder_dir, tumor_bam_rna, tumor_bai_rna]
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

    ch_versions = ch_versions.mix(NEO_ANNOTATE_FUSIONS.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, annotated_fusions ]
    ch_annotate_fusions_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(NEO_ANNOTATE_FUSIONS.out.annotated_fusions, ch_inputs),
            ch_isofox_inputs_sorted.skip.map { meta -> [meta, []] },
        )


    //
    // MODULE: Neo scorer
    //
    // Select input sources and prepare input channel
    // channel: [ meta_scorer, isofox_dir, purple_dir, sage_somatic_append, lilac_dir, neo_finder_dir, annotated_fusions ]
    ch_scorer_inputs = WorkflowOncoanalyser.groupByMeta(
        ch_isofox,
        ch_purple,
        ch_sage_somatic_append,
        ch_lilac,
        ch_finder_out,
        ch_annotate_fusions_out,
    )
        .map { meta, isofox_dir, purple_dir, sage_somatic_append, lilac_dir, neo_finder_dir, annotated_fusions ->

            def meta_scorer = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
                cancer_type: meta[Constants.InfoField.CANCER_TYPE],
            ]

            if (Utils.hasTumorRna(meta)) {
                meta_scorer.sample_rna_id = Utils.getTumorRnaSampleName(meta)
            }

            def inputs = [
                Utils.selectCurrentOrExisting(isofox_dir, meta, Constants.INPUT.ISOFOX_DIR),
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
                Utils.selectCurrentOrExisting(sage_somatic_append, meta, Constants.INPUT.SAGE_APPEND_VCF_TUMOR),
                Utils.selectCurrentOrExisting(lilac_dir, meta, Constants.INPUT.LILAC_DIR),
                neo_finder_dir,
                annotated_fusions,
            ]

            return [meta_scorer, *inputs]
        }

    // Run process
    NEO_SCORER(
        ch_scorer_inputs,
        ensembl_data_resources,
        neo_resources,
        cohort_tpm_medians,
    )

    ch_versions = ch_versions.mix(NEO_SCORER.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
