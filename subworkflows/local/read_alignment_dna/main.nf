//
// Align DNA reads
//

import Constants
import Utils

include { BWAMEM2_ALIGN  } from '../../../modules/local/bwa-mem2/mem/main'
include { FASTP          } from '../../../modules/local/fastp/main'

workflow READ_ALIGNMENT_DNA {
    take:
    // Sample data
    ch_inputs            // channel: [mandatory] [ meta ]

    // Reference data
    genome_fasta         // channel: [mandatory] /path/to/genome_fasta
    genome_bwamem2_index // channel: [mandatory] /path/to/genome_bwa-mem2_index_dir/

    // Params
    max_fastq_records    // numeric: [optional]  max number of FASTQ records per split
    umi_enable           // boolean: [mandatory] enable UMI processing
    umi_location         //  string: [optional]  fastp UMI location argument (--umi_loc)
    umi_length           // numeric: [optional]  fastp UMI length argument (--umi_len)
    umi_skip             // numeric: [optional]  fastp UMI skip argument (--umi_skip)

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Sort inputs, separate by tumor and normal
    // channel: [ meta ]
    ch_inputs_tumor_sorted = ch_inputs
        .branch { meta ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAM_DNA_TUMOR)
            runnable: Utils.hasTumorDnaFastq(meta) && !has_existing
            skip: true
        }

    ch_inputs_normal_sorted = ch_inputs
        .branch { meta ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAM_DNA_NORMAL)
            runnable: Utils.hasNormalDnaFastq(meta) && !has_existing
            skip: true
        }

    ch_inputs_donor_sorted = ch_inputs
        .branch { meta ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAM_DNA_DONOR)
            runnable: Utils.hasDonorDnaFastq(meta) && !has_existing
            skip: true
        }

    // Create FASTQ input channel
    // channel: [ meta_fastq, fastq_fwd, fastq_rev ]
    ch_fastq_inputs = Channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable.map { meta -> [meta, Utils.getTumorDnaSample(meta), 'tumor'] },
            ch_inputs_normal_sorted.runnable.map { meta -> [meta, Utils.getNormalDnaSample(meta), 'normal'] },
            ch_inputs_donor_sorted.runnable.map { meta -> [meta, Utils.getDonorDnaSample(meta), 'donor'] },
        )
        .flatMap { meta, meta_sample, sample_type ->
            meta_sample
                .getAt(Constants.FileType.FASTQ)
                .collect { key, fps ->
                    def (library_id, lane) = key

                    def meta_fastq = [
                        key: meta.group_id,
                        id: "${meta.group_id}_${meta_sample.sample_id}",
                        sample_id: meta_sample.sample_id,
                        library_id: library_id,
                        lane: lane,
                        sample_type: sample_type,
                    ]

                    return [meta_fastq, fps['fwd'], fps['rev']]
                }
        }

    //
    // MODULE: fastp
    //
    // Split FASTQ into chunks if requested for distributed processing
    // channel: [ meta_fastq_ready, fastq_fwd, fastq_fwd ]
    ch_fastqs_ready = Channel.empty()
    if (max_fastq_records > 0 || umi_enable) {

        // Run process
        FASTP(
            ch_fastq_inputs,
            max_fastq_records,
            umi_location,
            umi_length,
            umi_skip,
        )

        ch_versions = ch_versions.mix(FASTP.out.versions)

    }

    // Now prepare according to FASTQs splitting
    if (max_fastq_records > 0) {

        ch_fastqs_ready = FASTP.out.fastq
            .flatMap { meta_fastq, reads_fwd, reads_rev ->

                def data = [reads_fwd, reads_rev]
                    .transpose()
                    .collect { fwd, rev ->

                        def split_fwd = fwd.name.replaceAll('\\..+$', '')
                        def split_rev = rev.name.replaceAll('\\..+$', '')

                        assert split_fwd == split_rev

                        // NOTE(SW): split allows meta_fastq_ready to be unique, which is required during reunite below
                        def meta_fastq_ready = [
                            *:meta_fastq,
                            id: "${meta_fastq.id}_${split_fwd}",
                            split: split_fwd,
                        ]

                        return [meta_fastq_ready, fwd, rev]
                    }

                return data
            }

    } else {

        // Select appropriate source
        ch_fastq_source = umi_enable ? FASTP.out.fastq : ch_fastq_inputs

        ch_fastqs_ready = ch_fastq_source
            .map { meta_fastq, fastq_fwd, fastq_rev ->

                def meta_fastq_ready = [
                    *:meta_fastq,
                    split: null,
                ]

                return [meta_fastq_ready, fastq_fwd, fastq_rev]
            }

    }

    //
    // MODULE: BWA-MEM2
    //
    // Create process input channel
    // channel: [ meta_bwamem2, fastq_fwd, fastq_rev ]
    ch_bwamem2_inputs = ch_fastqs_ready
        .map { meta_fastq_ready, fastq_fwd, fastq_rev ->

            def meta_bwamem2 = [
                *:meta_fastq_ready,
                read_group: "${meta_fastq_ready.sample_id}.${meta_fastq_ready.library_id}.${meta_fastq_ready.lane}",
            ]

            return [meta_bwamem2, fastq_fwd, fastq_rev]
        }

    // Run process
    BWAMEM2_ALIGN(
        ch_bwamem2_inputs,
        genome_fasta,
        genome_bwamem2_index,
    )

    ch_versions = ch_versions.mix(BWAMEM2_ALIGN.out.versions)

    // Reunite BAMs
    // First, count expected BAMs per sample for non-blocking groupTuple op
    // channel: [ meta_count, group_size ]
    ch_sample_fastq_counts = ch_bwamem2_inputs
        .map { meta_bwamem2, reads_fwd, reads_rev ->

            def meta_count = [
                key: meta_bwamem2.key,
                sample_type: meta_bwamem2.sample_type,
            ]

            return [meta_count, meta_bwamem2]
        }
        .groupTuple()
        .map { meta_count, metas_bwamem2 -> return [meta_count, metas_bwamem2.size()] }

    // Now, group with expected size then sort into tumor and normal channels
    // channel: [ meta_group, [bam, ...], [bai, ...] ]
    ch_bams_united = ch_sample_fastq_counts
        .cross(
            // First element to match meta_count above for `cross`
            BWAMEM2_ALIGN.out.bam.map { meta_bwamem2, bam, bai -> [[key: meta_bwamem2.key, sample_type: meta_bwamem2.sample_type], bam, bai] }
        )
        .map { count_tuple, bam_tuple ->

            def group_size = count_tuple[1]
            def (meta_bam, bam, bai) = bam_tuple

            def meta_group = [
                *:meta_bam,
            ]

            return tuple(groupKey(meta_group, group_size), bam, bai)
        }
        .groupTuple()
        .branch { meta_group, bams, bais ->
            assert ['tumor', 'normal', 'donor'].contains(meta_group.sample_type)
            tumor: meta_group.sample_type == 'tumor'
            normal: meta_group.sample_type == 'normal'
            donor: meta_group.sample_type == 'donor'
            placeholder: true
        }

    // Set outputs, restoring original meta
    // channel: [ meta, [bam, ...], [bai, ...] ]
    ch_bam_tumor_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_bams_united.tumor, ch_inputs),
            ch_inputs_tumor_sorted.skip.map { meta -> [meta, [], []] },
        )

    ch_bam_normal_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_bams_united.normal, ch_inputs),
            ch_inputs_normal_sorted.skip.map { meta -> [meta, [], []] },
        )

    ch_bam_donor_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_bams_united.donor, ch_inputs),
            ch_inputs_donor_sorted.skip.map { meta -> [meta, [], []] },
        )

    emit:
    dna_tumor  = ch_bam_tumor_out  // channel: [ meta, [bam, ...], [bai, ...] ]
    dna_normal = ch_bam_normal_out // channel: [ meta, [bam, ...], [bai, ...] ]
    dna_donor  = ch_bam_donor_out  // channel: [ meta, [bam, ...], [bai, ...] ]

    versions   = ch_versions       // channel: [ versions.yml ]
}
