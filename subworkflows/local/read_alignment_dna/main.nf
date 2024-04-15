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
        ch_inputs              // channel: [mandatory] [ meta ]

        // Reference data
        genome_fasta           // channel: [mandatory] /path/to/genome_fasta
        genome_bwa_index       // channel: [mandatory] /path/to/genome_bwa_index_dir/
        genome_bwa_index_bseq  // channel: [mandatory] /path/to/genome_bwa_index_binary_seq
        genome_bwa_index_biidx // channel: [mandatory] /path/to/genome_bwa_index_bi-index

        // Params
        max_fastq_records      // numeric: [mandatory] max number of FASTQ records per split

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

        // Create FASTQ input channel
        // channel: [ meta_fastq, fastq_fwd, fastq_rev ]
        ch_fastq_inputs = Channel.empty()
            .mix(
                ch_inputs_tumor_sorted.runnable.map { meta -> [meta, Utils.getTumorDnaSample(meta), 'tumor'] },
                ch_inputs_normal_sorted.runnable.map { meta -> [meta, Utils.getNormalDnaSample(meta), 'normal'] },
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
        if (max_fastq_records > 0) {

            // Run process
            FASTP(
                ch_fastq_inputs,
                max_fastq_records,
            )

            ch_versions = ch_versions.mix(FASTP.out.versions)

            // Prepare outputs within conditional block
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

            ch_fastqs_ready = ch_fastq_inputs
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
        // channel: [ meta_bwa, fastq_fwd, fastq_rev ]
        ch_bwa_inputs = ch_fastqs_ready
            .map { meta_fastq_ready, fastq_fwd, fastq_rev ->

                def meta_bwa = [
                    *:meta_fastq_ready,


                    // TODO(SW): understand target format
                    read_group: "${meta_fastq_ready.sample_id}.${meta_fastq_ready.library_id}.${meta_fastq_ready.lane}",


                ]

                return [meta_bwa, fastq_fwd, fastq_rev]
            }

        // Run process
        BWAMEM2_ALIGN(
            ch_bwa_inputs,
            genome_fasta,
            genome_bwa_index,
            genome_bwa_index_bseq,
            genome_bwa_index_biidx,
        )

        ch_versions = ch_versions.mix(BWAMEM2_ALIGN.out.versions)

        // Reunite BAMs
        // First, count expected BAMs per sample for non-blocking groupTuple op
        // channel: [ meta_count, group_size ]
        ch_sample_fastq_counts = ch_bwa_inputs
            .map { meta_bwa, reads_fwd, reads_rev ->

                def meta_count = [
                    key: meta_bwa.key,
                    sample_type: meta_bwa.sample_type,
                ]

                return [meta_count, meta_bwa]
            }
            .groupTuple()
            .map { meta_count, meta_bwas -> return [meta_count, meta_bwas.size()] }

        // Now, group with expected size then sort into tumor and normal channels
        // channel: [ meta_group, [bam, ...], [bai, ...] ]
        ch_bams_united = ch_sample_fastq_counts
            .cross(
                // First element to match meta_count above for `cross`
                BWAMEM2_ALIGN.out.bam.map { meta_bwa, bam, bai -> [[key: meta_bwa.key, sample_type: meta_bwa.sample_type], bam, bai] }
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
                assert ['tumor', 'normal'].contains(meta_group.sample_type)
                tumor: meta_group.sample_type == 'tumor'
                normal: meta_group.sample_type == 'normal'
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

    emit:
        dna_tumor  = ch_bam_tumor_out  // channel: [ meta, [bam, ...], [bai, ...] ]
        dna_normal = ch_bam_normal_out // channel: [ meta, [bam, ...], [bai, ...] ]
        versions   = ch_versions       // channel: [ versions.yml ]
}
