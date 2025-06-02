//
// Bam Tools calculates summary statistics for BAMs
//

import Constants
import Utils

include { BAMTOOLS } from '../../../modules/local/bamtools/main'

workflow BAMTOOLS_METRICS {
    take:
    // Sample data
    ch_inputs      // channel: [mandatory] [ meta ]
    ch_tumor_bam   // channel: [mandatory] [ meta, bam, bai ]
    ch_normal_bam  // channel: [mandatory] [ meta, bam, bai ]

    // Reference data
    genome_fasta   // channel: [mandatory] /path/to/genome_fasta
    genome_version // channel: [mandatory] genome version

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Sort inputs, separate by tumor and normal
    // channel: runnable: [ meta, bam, bai ]
    // channel: skip: [ meta ]
    ch_inputs_tumor_sorted = ch_tumor_bam
        .map { meta, bam, bai ->
            return [
                meta,
                Utils.selectCurrentOrExisting(bam, meta, Constants.INPUT.BAM_REDUX_DNA_TUMOR),
                bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_TUMOR),
            ]
        }
        .branch { meta, bam, bai ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAMTOOLS_DIR_TUMOR)
            runnable: bam && !has_existing
            skip: true
                return meta
        }

    // channel: runnable: [ meta, bam, bai ]
    // channel: skip: [ meta ]
    ch_inputs_normal_sorted = ch_normal_bam
        .map { meta, bam, bai ->
            return [
                meta,
                Utils.selectCurrentOrExisting(bam, meta, Constants.INPUT.BAM_REDUX_DNA_NORMAL),
                bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_NORMAL),
            ]
        }
        .branch { meta, bam, bai ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAMTOOLS_DIR_NORMAL)
            runnable: bam && !has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_bamtools, bam, bai ]
    ch_bamtools_inputs = Channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable.map { meta, bam, bai -> [meta, Utils.getTumorDnaSample(meta), 'tumor', bam, bai] },
            ch_inputs_normal_sorted.runnable.map { meta, bam, bai -> [meta, Utils.getNormalDnaSample(meta), 'normal', bam, bai] },
        )
        .map { meta, meta_sample, sample_type, bam, bai ->

            def meta_bamtools = [
                key: meta.group_id,
                id: "${meta.group_id}_${meta_sample.sample_id}",
                sample_id: meta_sample.sample_id,
                sample_type: sample_type,
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

    // Sort into a tumor and normal channel
    ch_bamtools_out = BAMTOOLS.out.metrics_dir
        .branch { meta_bamtools, metrics_dir ->
            assert ['tumor', 'normal'].contains(meta_bamtools.sample_type)
            tumor: meta_bamtools.sample_type == 'tumor'
            normal: meta_bamtools.sample_type == 'normal'
            placeholder: true
        }

    // Set outputs, restoring original meta
    // channel: [ meta, metrics_dir ]
    ch_somatic_metrics_dir = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_bamtools_out.tumor, ch_inputs),
            ch_inputs_tumor_sorted.skip.map { meta -> [meta, []] },
        )

    ch_germline_metrics_dir = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_bamtools_out.normal, ch_inputs),
            ch_inputs_normal_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    somatic  = ch_somatic_metrics_dir  // channel: [ meta, metrics_dir ]
    germline = ch_germline_metrics_dir // channel: [ meta, metrics_dir ]

    versions = ch_versions             // channel: [ versions.yml ]
}
