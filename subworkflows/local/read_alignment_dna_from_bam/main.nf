//
// Align DNA reads from an existing BAM/CRAM by streaming through samtools fastq
// Note: this does not include the option to chunk up the bam, as is available in the FASTQ-based workflow.
//

import Constants
import Utils

include { BWAMEM2_ALIGN_FROM_BAM } from '../../../modules/local/bwa-mem2/mem_from_bam/main'

workflow READ_ALIGNMENT_DNA_FROM_BAM {
    take:
    // Sample data
    ch_inputs            // channel: [mandatory] [ meta ]

    // Reference data
    genome_fasta         // channel: [mandatory] /path/to/genome_fasta
    genome_bwamem2_index // channel: [mandatory] /path/to/genome_bwa-mem2_index_dir/

    main:
    ch_versions = Channel.empty()

    // Branch inputs: process samples that have BAM/CRAM inputs
    ch_inputs_tumor_sorted = ch_inputs
        .branch { meta ->
            runnable: Utils.hasTumorDnaBam(meta)
            skip: true
        }

    ch_inputs_normal_sorted = ch_inputs
        .branch { meta ->
            runnable: Utils.hasNormalDnaBam(meta)
            skip: true
        }

    ch_inputs_donor_sorted = ch_inputs
        .branch { meta ->
            runnable: Utils.hasDonorDnaBam(meta)
            skip: true
        }

    // Build a flat channel of [ meta_bwamem2, bam ] from input BAM/CRAM.
    ch_bwamem2_inputs = Channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable.map { meta ->
                def meta_sample = Utils.getTumorDnaSample(meta)
                def sample_id   = meta_sample.getOrDefault('longitudinal_sample_id', meta_sample['sample_id'])
                def meta_bwamem2 = [
                    key:         meta.group_id,
                    id:          "${meta.group_id}_${sample_id}",
                    sample_id:   sample_id,
                    sample_type: 'tumor',
                    read_group:  "${sample_id}.realign",
                ]
                def bam = Utils.getTumorDnaBam(meta)
                assert bam != null && bam != [] : "Missing tumor DNA BAM/CRAM for ${meta.group_id}"
                return [meta_bwamem2, bam]
            },
            ch_inputs_normal_sorted.runnable.map { meta ->
                def meta_sample = Utils.getNormalDnaSample(meta)
                def sample_id   = meta_sample['sample_id']
                def meta_bwamem2 = [
                    key:         meta.group_id,
                    id:          "${meta.group_id}_${sample_id}",
                    sample_id:   sample_id,
                    sample_type: 'normal',
                    read_group:  "${sample_id}.realign",
                ]
                def bam = Utils.getNormalDnaBam(meta)
                assert bam != null && bam != [] : "Missing normal DNA BAM/CRAM for ${meta.group_id}"
                return [meta_bwamem2, bam]
            },
            ch_inputs_donor_sorted.runnable.map { meta ->
                def meta_sample = Utils.getDonorDnaSample(meta)
                def sample_id   = meta_sample['sample_id']
                def meta_bwamem2 = [
                    key:         meta.group_id,
                    id:          "${meta.group_id}_${sample_id}",
                    sample_id:   sample_id,
                    sample_type: 'donor',
                    read_group:  "${sample_id}.realign",
                ]
                def bam = Utils.getDonorDnaBam(meta)
                assert bam != null && bam != [] : "Missing donor DNA BAM/CRAM for ${meta.group_id}"
                return [meta_bwamem2, bam]
            },
        )
        .filter { meta_bwamem2, bam -> bam != null && bam != [] }

    BWAMEM2_ALIGN_FROM_BAM(
        ch_bwamem2_inputs,
        genome_fasta,
        genome_bwamem2_index,
    )

    ch_versions = ch_versions.mix(BWAMEM2_ALIGN_FROM_BAM.out.versions)

    // Reunite BAMs (count = 1 per sample for simple single-BAM-in case)
    ch_sample_counts = ch_bwamem2_inputs
        .map { meta_bwamem2, bam -> [[key: meta_bwamem2.key, sample_type: meta_bwamem2.sample_type], meta_bwamem2] }
        .groupTuple()
        .map { meta_count, metas -> [meta_count, metas.size()] }

    ch_bams_united = ch_sample_counts
        .cross(
            BWAMEM2_ALIGN_FROM_BAM.out.bam
                .map { meta_bwamem2, bam, bai -> [[key: meta_bwamem2.key, sample_type: meta_bwamem2.sample_type], bam, bai] }
        )
        .map { count_tuple, bam_tuple ->
            def group_size = count_tuple[1]
            def (meta_bam, bam, bai) = bam_tuple
            return tuple(groupKey([*:meta_bam], group_size), bam, bai)
        }
        .groupTuple()
        .branch { meta_group, bams, bais ->
            assert ['tumor', 'normal', 'donor'].contains(meta_group.sample_type)
            tumor:  meta_group.sample_type == 'tumor'
            normal: meta_group.sample_type == 'normal'
            donor:  meta_group.sample_type == 'donor'
            placeholder: true
        }

    ch_bam_tumor_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_bams_united.tumor,  ch_inputs),
            ch_inputs_tumor_sorted.skip.map  { meta -> [meta, [], []] },
        )

    ch_bam_normal_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_bams_united.normal, ch_inputs),
            ch_inputs_normal_sorted.skip.map { meta -> [meta, [], []] },
        )

    ch_bam_donor_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_bams_united.donor,  ch_inputs),
            ch_inputs_donor_sorted.skip.map  { meta -> [meta, [], []] },
        )

    emit:
    dna_tumor  = ch_bam_tumor_out   // channel: [ meta, [bam, ...], [bai, ...] ]
    dna_normal = ch_bam_normal_out  // channel: [ meta, [bam, ...], [bai, ...] ]
    dna_donor  = ch_bam_donor_out   // channel: [ meta, [bam, ...], [bai, ...] ]

    versions   = ch_versions        // channel: [ versions.yml ]
}
