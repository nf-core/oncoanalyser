//
// XXX
//

import Constants
import Utils

include { ISOFOX               } from '../../modules/local/isofox/main'
include { LILAC                } from '../../modules/local/lilac/main'
include { NEO_FINDER           } from '../../modules/local/neo/finder/main'
include { NEO_SCORER           } from '../../modules/local/neo/scorer/main'

workflow NEO_PREDICTION{
    take:
        // Sample data
        ch_inputs              // channel: [mandatory] [ meta ]
        ch_isofox              // channel: [mandatory] [ meta, isofox_dir ]
        ch_purple              // channel: [mandatory] [ meta, purple_dir ]
        ch_sage_somatic_append // channel: [mandatory] [ meta, sage_append_vcf ]
        ch_lilac               // channel: [mandatory] [ meta, lilac_dir ]
        ch_linx                // channel: [mandatory] [ meta, linx_dir ]

        // Reference data
        genome_version         // channel: [mandatory] genome version
        genome_fasta           // channel: [mandatory] /path/to/genome_fasta
        genome_fai             // channel: [mandatory] /path/to/genome_fai
        ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/
        neo_resources          // channel: [mandatory] /path/to/neo_resources/
        cohort_tpm_medians     // channel: [mandatory] /path/to/cohort_tpm_medians/

    main:
        // Channel for versions.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Step 1: Identify neoepitopes from Purple somatic variants and Linx's (neoepitope) fusions

        // Select input sources
        // channel: [ meta, isofox_dir, purple_dir, linx_annotation_dir ]
        ch_inputs_finder_selected = WorkflowOncoanalyser.groupByMeta(
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
        ch_inputs_finder_sorted = ch_inputs_finder_selected
            .branch { meta, purple_dir, linx_annotation_dir ->

                def has_normal_dna = Utils.hasNormalDnaBam(meta)

                def has_runnable_inputs = purple_dir && linx_annotation_dir && has_normal_dna

                runnable: has_runnable_inputs
                skip: true
                    return meta
            }

        // Create process input channel
        // channel: sample_data: [ meta, purple_dir, linx_annotation_dir ]
        ch_finder_inputs = ch_inputs_finder_sorted.runnable
            .map{ meta, purple_dir, linx_annotation_dir ->

                def meta_neo_finder = [
                    key: meta.group_id,
                    id: meta.group_id,
                    sample_id: Utils.getTumorDnaSampleName(meta),
                ]

                return [meta_neo_finder, purple_dir, linx_annotation_dir]
            }


        // Feeding the Neo process raw inputs for demo purposes only
        NEO_FINDER(
            ch_finder_inputs,
            genome_fasta,
            genome_version,
            ensembl_data_resources,
        )

        ch_versions = ch_versions.mix(NEO_FINDER.out.versions)

        // Set outputs, restoring original meta
        // channel: [ meta, neo_finder_dir ]
        ch_finder_outputs = WorkflowOncoanalyser.restoreMeta(NEO_FINDER.out.neo_finder_dir, ch_inputs)

        // Step 2: When RNA is present, annotate the fusion-derived neoepitope with RNA using Isofox

        /*

        // Select input sources
        // channel: [ meta, neo_finder_dir, tumor_bam_rna, tumor_bai_rna ]
        ch_inputs_isofox_sorted = WorkflowOncoanalyser.groupByMeta(
            ch_finder_outputs,
            // channel: [ meta, tumor_rna_bam (optional), tumor_rna_bai (optional) ]
            ch_inputs
                .map { meta ->
                    def has_rna = Utils.hasTumorRnaBam(meta)

                    return [
                        meta,
                        has_rna ? Utils.getTumorRnaBam(meta) : [],
                        has_rna ? Utils.getTumorRnaBai(meta) : [],
                    ]
                },

            )

        // Sort inputs
        ch_inputs_isofox_sorted = ch_finder_outputs
            .branch {

                def has_rna = Utils.hasTumorRnaBam(meta)



                runnable:
                skip:
                    meta

            }

        // Create process input channel
        // channel: [ meta_isofox, neo_finder_dir, tumor_bam_rna, tumor_bai_rna ]
        ch_isofox_inputs = ch_inputs_isofox_sorted.runnable
            .map { meta, neo_finder_dir, tumor_bam_rna, tumor_bai_rna ->

                def meta_isofox = [
                    key: meta.group_id,
                    id: meta.group_id,
                    sample_id: Utils.getTumorDnaSampleName(meta),
                ]

                return [meta_isofox, Utils.getTumorRnaBam(meta), Utils.getTumorRnaBai(meta)]
            }

        // Run process
        ISOFOX_NEO(
            ch_isofox_inputs,
            isofox_read_length,
            genome_fasta,
            genome_version,
            genome_fai,
            ensembl_data_resources,
        )

        ch_versions = ch_versions.mix(ISOFOX.out.versions)

        // Set outputs, restoring original meta
        // channel: [ meta, isofox_dir ]
        ch_outputs = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(ISOFOX.out.isofox_neo_dir, ch_inputs),
                ch_inputs_sorted.skip.map { meta -> [meta, []] },

        */

        // ch_finder_outputs

        // Step 3: Run Neo's binding prediction routine for neoepitope's pHLAs, taking in Lilac HLA alleles and previously
        // derived neoepitopes with RNA annotation if it was available

        // Select input sources
        // channel: [ meta, isofox_dir, purple_dir, lilac_dir, isofox_dir ]
        // TO_DO - how to pass in the directories from step 1 and 2 (if run) above
        ch_inputs_scorer_selected = WorkflowOncoanalyser.groupByMeta(
            ch_purple,
            ch_linx,
            ch_isofox,
        )
            .map { meta, purple_dir, lilac_dir ->

                def inputs = [
                    Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
                    Utils.selectCurrentOrExisting(lilac_dir, meta, Constants.INPUT.LILAC),
                    Utils.selectCurrentOrExisting(isofox_dir, meta, Constants.INPUT.ISOFOX),
                ]

                return [meta, *inputs]
            }

        // Sort inputs
        // channel: runnable: [ meta, purple_dir, lilac_dir,isofox_dir ]
        // channel: skip: [ meta ]
        ch_inputs_scorer_sorted = ch_inputs_scorer_selected
            .branch { meta, purple_dir, lilac_dir, isofox_dir ->

                def has_normal_dna = Utils.hasNormalDnaBam(meta)

                def has_runnable_inputs = purple_dir && lilac_dir && has_normal_dna

                runnable: has_runnable_inputs
                skip: true
                    return meta
            }

        // Create process input channel
        // channel: sample_data: [ meta, purple_dir, linx_annotation_dir ]
        ch_scorer_inputs = ch_inputs_scorer_sorted.runnable
            .map{ meta, purple_dir, linx_annotation_dir ->

                def meta_neo_scorer = [
                    key: meta.group_id,
                    id: meta.group_id,
                    sample_id: Utils.getTumorDnaSampleName(meta),
                ]

                return [meta_neo_scorer, purple_dir, lilac_dir, isofox_dir]
            }


        // Feeding the Neo process raw inputs for demo purposes only
        NEO_SCORER(
            ch_scorer_inputs,
            genome_fasta,
            genome_version,
            ensembl_data_resources,
            neo_resources,
            cohort_tpm_medians
        )

        ch_versions = ch_versions.mix(NEO_SCORER.out.versions)

        // Set outputs, restoring original meta
        // channel: [ meta, neo_scorer_dir ]
        ch_scorer_outputs = WorkflowOncoanalyser.restoreMeta(NEO_SCORER.out.neo_scorer_dir, ch_inputs)

    emit:
        versions = ch_versions // channel: [ versions.yml ]
}
