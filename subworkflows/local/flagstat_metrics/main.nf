//
// SAMtools flagstat generates statistics for read alignments from the SAM FLAG field
//

import Constants
import Utils

include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'

workflow FLAGSTAT_METRICS {
    take:
        // Sample data
        ch_inputs     // channel: [mandatory] [ meta ]
        ch_tumor_bam  // channel: [mandatory] [ meta, bam, bai ]
        ch_normal_bam // channel: [mandatory] [ meta, bam, bai ]

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
                    Utils.selectCurrentOrExisting(bam, meta, Constants.INPUT.BAM_MARKDUPS_DNA_TUMOR),
                    bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_TUMOR),
                ]
            }
            .branch { meta, bam, bai ->
                def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.FLAGSTAT_TUMOR)
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
                    Utils.selectCurrentOrExisting(bam, meta, Constants.INPUT.BAM_MARKDUPS_DNA_NORMAL),
                    bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_NORMAL),
                ]
            }
            .branch { meta, bam, bai ->
                def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.FLAGSTAT_NORMAL)
                runnable: bam && !has_existing
                skip: true
                    return meta
            }

        // Create process input channel
        // channel: [ meta_flagstat, bam, bai ]
        ch_flagstat_inputs = Channel.empty()
            .mix(
                ch_inputs_tumor_sorted.runnable.map { meta, bam, bai -> [meta, Utils.getTumorDnaSample(meta), 'tumor', bam, bai] },
                ch_inputs_normal_sorted.runnable.map { meta, bam, bai -> [meta, Utils.getNormalDnaSample(meta), 'normal', bam, bai] },
            )
            .map { meta, meta_sample, sample_type, bam, bai ->

                def meta_flagstat = [
                    key: meta.group_id,
                    id: "${meta.group_id}_${meta_sample.sample_id}",
                    sample_id: meta_sample.sample_id,
                    sample_type: sample_type,
                ]

                return [meta_flagstat, bam, bai]
            }

        // Run process
        SAMTOOLS_FLAGSTAT(
            ch_flagstat_inputs,
        )

        ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

        // Sort into a tumor and normal channel
        ch_flagstat_out = SAMTOOLS_FLAGSTAT.out.flagstat
            .branch { meta_flagstat, flagstat ->
                assert ['tumor', 'normal'].contains(meta_flagstat.sample_type)
                tumor: meta_flagstat.sample_type == 'tumor'
                normal: meta_flagstat.sample_type == 'normal'
                placeholder: true
            }

        // Set outputs, restoring original meta
        // channel: [ meta, flagstat ]
        ch_somatic_flagstat = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(ch_flagstat_out.tumor, ch_inputs),
                ch_inputs_tumor_sorted.skip.map { meta -> [meta, []] },
            )

        ch_germline_flagstat = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(ch_flagstat_out.normal, ch_inputs),
                ch_inputs_normal_sorted.skip.map { meta -> [meta, []] },
            )

    emit:
        somatic  = ch_somatic_flagstat  // channel: [ meta, flagstat ]
        germline = ch_germline_flagstat // channel: [ meta, flagstat ]

        versions = ch_versions          // channel: [ versions.yml ]
}
