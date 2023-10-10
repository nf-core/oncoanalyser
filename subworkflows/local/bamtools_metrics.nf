//
// Bam Tools calculates summary statistics for BAMs
//

import Constants
import Utils

include { BAMTOOLS } from '../../modules/local/bamtools/main'

workflow BAMTOOLS_METRICS {
    take:
        // Sample data
        ch_inputs      // channel: [mandatory] [ meta ]

        // Reference data
        genome_fasta   // channel: [mandatory] /path/to/genome_fasta
        genome_version // channel: [mandatory] genome version

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Sort inputs
        // channel: [ meta ]
        ch_inputs_sorted = ch_inputs
            .branch { meta ->

                def has_tumor_dna = Utils.hasTumorDnaBam(meta)
                def has_normal_dna = Utils.hasNormalDnaBam(meta)

                runnable: has_tumor_dna || has_normal_dna
                skip: true
            }

        // Flatten into BAM/BAI pairs, select inputs that are eligible to run
        // channel: runnable: [ meta_extra, bam, bai ]
        // channel: skip: [ meta_extra ]
        ch_bams_bais_sorted = ch_inputs_sorted.runnable
            .flatMap { meta ->

                def tumor_sample_id = []
                def tumor_bam = []
                def tumor_bai = []

                def normal_sample_id = []
                def normal_bam = []
                def normal_bai = []


                if (Utils.hasTumorDnaBam(meta)) {
                    tumor_sample_id = Utils.getTumorDnaSampleName(meta)
                    tumor_bam = Utils.getTumorDnaBam(meta)
                    tumor_bai = Utils.getTumorDnaBai(meta)
                }

                if (Utils.hasNormalDnaBam(meta)) {
                    normal_sample_id = Utils.getNormalDnaSampleName(meta)
                    normal_bam = Utils.getNormalDnaBam(meta)
                    normal_bai = Utils.getNormalDnaBai(meta)
                }

                return [
                    [[key: meta.group_id, *:meta, sample_id: tumor_sample_id, sample_type: 'tumor'], tumor_bam, tumor_bai],
                    [[key: meta.group_id, *:meta, sample_id: normal_sample_id, sample_type: 'normal'], normal_bam, normal_bai],
                ]
            }
            .branch { meta_extra, bam, bai ->

                def input_key
                if (meta_extra.sample_type == 'tumor') {
                    input_key = Constants.INPUT.BAMTOOLS_TUMOR
                } else if (meta_extra.sample_type == 'normal') {
                    input_key = Constants.INPUT.BAMTOOLS_NORMAL
                } else {
                    assert false
                }

                def has_existing = Utils.hasExistingInput(meta_extra, input_key)

                runnable: bam && bai && !has_existing
                skip: true
                    return meta_extra
            }

        // Create process input channel
        // channel: [ meta_bamtools, bam, bai ]
        ch_bamtools_inputs = ch_bams_bais_sorted.runnable
            .map { meta_extra, bam, bai ->

                def meta_bamtools = [
                    key: meta_extra.group_id,
                    id: "${meta_extra.group_id}__${meta_extra.sample_id}",
                    sample_type: meta_extra.sample_type,
                ]

                return [meta_bamtools, bam, bai]
            }

        // Run process
        BAMTOOLS(
            ch_bamtools_inputs,
            genome_fasta,
            genome_version,
        )

        ch_versions = ch_versions.mix(BAMTOOLS.out.versions)

        // Sort outputs into tumor and normal channels, adding partial skip entries
        // channel: [ meta_bamtools, metrics ]
        ch_outputs_sorted = Channel.empty()
            .mix(
                BAMTOOLS.out.metrics,
                ch_bams_bais_sorted.skip.map { meta -> [meta, []] },
            )
            .branch { meta_bamtools, metrics ->
                tumor: meta_bamtools.sample_type == 'tumor'
                normal: meta_bamtools.sample_type == 'normal'
            }

        // Set outputs, restoring original meta, including full skip entries
        ch_somatic_metrics = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(ch_outputs_sorted.tumor, ch_inputs),
                ch_inputs_sorted.skip.map { meta -> [meta, []] },
            )

        ch_germline_metrics = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(ch_outputs_sorted.normal, ch_inputs),
                ch_inputs_sorted.skip.map { meta -> [meta, []] },
            )

    emit:
        somatic  = ch_somatic_metrics  // channel: [ meta, metrics ]
        germline = ch_germline_metrics // channel: [ meta, metrics ]

        versions = ch_versions         // channel: [ versions.yml ]
}
