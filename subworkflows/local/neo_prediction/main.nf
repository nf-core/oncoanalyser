//
// Neo identifies and scores neoepitopes
//

include { NEO_ANNOTATE_FUSIONS } from '../../../modules/local/neo/annotate_fusions/main'
include { NEO_FINDER           } from '../../../modules/local/neo/finder/main'
include { NEO_SCORER           } from '../../../modules/local/neo/scorer/main'

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
    ch_finder_inputs_selected = WorkflowOncoanalyser.groupByMeta(
        ch_purple_dir,
        ch_linx_annotation_dir,
    )
        .map { meta, purple_dir, linx_annotation_dir ->

            return [
                meta,
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
                Utils.selectCurrentOrExisting(linx_annotation_dir, meta, Constants.INPUT.LINX_ANNO_DIR_TUMOR),
            ]

        }

    // Sort inputs
    // channel: runnable: [ meta, purple_dir, linx_annotation_dir ]
    // channel: skip: [ meta ]
    ch_finder_inputs_sorted = ch_finder_inputs_selected
        .branch { meta, purple_dir, linx_annotation_dir ->

            def has_normal_dna = Utils.hasNormalDna(meta)

            def tumor_id = Utils.getTumorDnaSampleName(meta)
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

    // Set outputs, restoring original meta
    // channel: [ meta, neo_finder_dir ]
    ch_finder_out = WorkflowOncoanalyser.restoreMeta(channel.topic('neo_finder_dir'), ch_inputs)

    //
    // MODULE: Fusion annotation (Isofox)
    //
    // Annotate the fusion-derived neoepitope using Isofox where RNA data is available

    // Select input sources and sort
    // channel: runnable: [ meta, neo_finder_dir, tumor_rna_aln, tumor_rna_idx ]
    // channel: skip: [ meta ]
    ch_isofox_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_finder_out,
        ch_tumor_rna_aln,
    )
        .map { meta, neo_finder_dir, tumor_rna_aln, tumor_rna_idx ->
            return [
                meta,
                neo_finder_dir,
                Utils.selectCurrentOrExisting(tumor_rna_aln, meta, Constants.INPUT.ALN_RNA_TUMOR),
                Utils.selectCurrentOrExisting(tumor_rna_idx, meta, Constants.INPUT.IDX_RNA_TUMOR),
            ]
        }
        .branch { meta, neo_finder_dir, tumor_rna_aln, tumor_rna_idx ->
            runnable: Utils.hasTumorRna(meta)
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_isofox, neo_finder_dir, tumor_rna_aln, tumor_rna_idx ]
    ch_isofox_inputs = ch_isofox_inputs_sorted.runnable
        .map { meta, neo_finder_dir, tumor_rna_aln, tumor_rna_idx ->

            def meta_isofox = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
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
            WorkflowOncoanalyser.restoreMeta(channel.topic('neo_annotated_fusions_tsv'), ch_inputs),
            ch_isofox_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    //
    // MODULE: Neo scorer
    //
    // Select input sources and prepare input channel
    // channel: [ meta_scorer, isofox_dir, purple_dir, sage_append_dir_somatic, lilac_dir, neo_finder_dir, annotated_fusions ]
    ch_scorer_inputs = WorkflowOncoanalyser.groupByMeta(
        ch_isofox_dir,
        ch_purple_dir,
        ch_sage_append_dir_somatic,
        ch_lilac_dir,
        ch_finder_out,
        ch_annotate_fusions_out,
    )
        .map { meta, isofox_dir, purple_dir, sage_append_dir_somatic, lilac_dir, neo_finder_dir, annotated_fusions ->

            def meta_scorer = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
                cancer_type: meta[Constants.InfoField.CANCER_TYPE],
            ]

            def sage_append_vcf_somatic = []
            if (Utils.hasTumorRna(meta)) {
                meta_scorer.sample_rna_id = Utils.getTumorRnaSampleName(meta)

                def sage_append_dir_somatic_selected = Utils.selectCurrentOrExisting(sage_append_dir_somatic, meta, Constants.INPUT.SAGE_APPEND_DIR_TUMOR)
                sage_append_vcf_somatic = file(sage_append_dir_somatic_selected).resolve("${meta_scorer.sample_id}.sage.append.vcf.gz")
            }

            return [
                meta_scorer,
                Utils.selectCurrentOrExisting(isofox_dir, meta, Constants.INPUT.ISOFOX_DIR),
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
                sage_append_vcf_somatic,
                Utils.selectCurrentOrExisting(lilac_dir, meta, Constants.INPUT.LILAC_DIR),
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
}
