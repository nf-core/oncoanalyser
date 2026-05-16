//
// NEO identifies and scores neoepitopes
//

import Constants
import Utils

include { NEO_FINDER            } from '../../../modules/local/neo/finder/main'
include { NEO_ANNOTATE_FUSIONS  } from '../../../modules/local/neo/annotate_fusions/main'
include { NEO_SCORER            } from '../../../modules/local/neo/scorer/main'

workflow NEO_PREDICTION {
    take:
    // Sample data
    ch_inputs             // channel: [mandatory] [ meta ]
    ch_tumor_rna_bam      // channel: [mandatory] [ meta, bam, bai ]
    ch_isofox             // channel: [mandatory] [ meta, isofox_dir ]
    ch_purple             // channel: [mandatory] [ meta, purple_dir ]
    ch_sage_somatic_append // channel: [mandatory] [ meta, sage_append_dir ]
    ch_lilac              // channel: [mandatory] [ meta, lilac_dir ]
    ch_linx_somatic       // channel: [mandatory] [ meta, linx_somatic_dir ]

    // Reference data
    genome_fasta           // channel: [mandatory] /path/to/genome_fasta
    genome_version         // channel: [mandatory] genome version
    genome_fai             // channel: [mandatory] /path/to/genome_fai
    ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/
    neo_resources          // channel: [mandatory] /path/to/neo_resources/
    cohort_tpm_medians     // channel: [mandatory] /path/to/cohort_tpm_medians

    // Params
    isofox_read_length     // integer: [mandatory] Isofox read length

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    //
    // MODULE: Neo finder
    //
    // Select input sources and sort
    // channel: runnable: [ meta, purple_dir, lilac_dir, linx_somatic_dir ]
    // channel: skip: [ meta ]
    ch_finder_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_purple,
        ch_lilac,
        ch_linx_somatic,
    )
        .map { meta, purple_dir, lilac_dir, linx_somatic_dir ->
            return [
                meta,
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
                Utils.selectCurrentOrExisting(lilac_dir, meta, Constants.INPUT.LILAC_DIR),
                Utils.selectCurrentOrExisting(linx_somatic_dir, meta, Constants.INPUT.LINX_ANNO_DIR_TUMOR),
            ]
        }
        .branch { meta, purple_dir, lilac_dir, linx_somatic_dir ->
            runnable: Utils.hasTumorDna(meta)
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_finder, purple_dir, lilac_dir, linx_somatic_dir ]
    ch_finder_inputs = ch_finder_inputs_sorted.runnable
        .map { meta, purple_dir, lilac_dir, linx_somatic_dir ->

            def meta_finder = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta, primary: true),
            ]

            return [meta_finder, purple_dir, lilac_dir, linx_somatic_dir]
        }

    // Run process
    NEO_FINDER(
        ch_finder_inputs,
        neo_resources,
    )

    ch_versions = ch_versions.mix(NEO_FINDER.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, neo_finder_dir ]
    ch_finder_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(NEO_FINDER.out.neo_finder_dir, ch_inputs),
            ch_finder_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    //
    // MODULE: Annotate fusions
    //
    // Select input sources and sort for fusion annotation
    // channel: runnable: [ meta, neo_finder_dir, tumor_bam, tumor_bai ]
    // channel: skip: [ meta ]
    ch_isofox_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_finder_out,
        ch_tumor_rna_bam,
    )
        .map { meta, neo_finder_dir, tumor_bam, tumor_bai ->
            return [
                meta,
                neo_finder_dir,
                params.realign_bam
                    ? tumor_bam
                    : Utils.selectCurrentOrExisting(tumor_bam, meta, Constants.INPUT.BAM_RNA_TUMOR),
                params.realign_bam
                    ? tumor_bai
                    : Utils.selectCurrentOrExisting(tumor_bai, meta, Constants.INPUT.BAI_RNA_TUMOR),
            ]
        }
        .branch { meta, neo_finder_dir, tumor_bam, tumor_bai ->
            runnable: Utils.hasTumorRna(meta)
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
                sample_id: Utils.getTumorDnaSampleName(meta, primary: true),
                cancer_type: meta[Constants.InfoField.CANCER_TYPE],
            ]

            def sage_somatic_append_vcf = []
            if (Utils.hasTumorRna(meta)) {
                meta_scorer.sample_rna_id = Utils.getTumorRnaSampleName(meta)

                def sage_somatic_append_selected = Utils.selectCurrentOrExisting(sage_somatic_append, meta, Constants.INPUT.SAGE_APPEND_DIR_TUMOR)
                sage_somatic_append_vcf = file(sage_somatic_append_selected).resolve("${meta_scorer.sample_id}.sage.append.vcf.gz")
            }

            def inputs = [
                Utils.selectCurrentOrExisting(isofox_dir, meta, Constants.INPUT.ISOFOX_DIR),
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
                sage_somatic_append_vcf,
                Utils.selectCurrentOrExisting(lilac_dir, meta, Constants.INPUT.LILAC_DIR),
                neo_finder_dir,
                annotated_fusions,
            ]

            return [meta_scorer, *inputs]
        }

    // Run process
    NEO_SCORER(
        ch_scorer_inputs,
        neo_resources,
        cohort_tpm_medians,
    )

    ch_versions = ch_versions.mix(NEO_SCORER.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
