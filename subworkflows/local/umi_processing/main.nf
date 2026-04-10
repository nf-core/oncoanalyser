//
// Align DNA reads
//

include { FASTQ_TOOLS } from '../../../modules/local/fastq-tools/main'

workflow UMI_PROCESSING {
    take:
    // Sample data
    ch_inputs            // channel: [mandatory] [ meta ]

    // Reference data
    known_umis           // channel: [mandatory] /path/to/known_umis_file/

    // Params
    umi_location         // string:  [optional]  fastp UMI location argument (--umi_loc)
    umi_length           // numeric: [optional]  fastp UMI length argument (--umi_len)
    umi_skip             // numeric: [optional]  fastp UMI skip argument (--umi_skip)

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    ch_inputs.map { meta ->

        def meta_sample_entries = [
            tumor: Inputs.getTumorDnaSample(meta),
            normal: Inputs.getNormalDnaSample(meta),
            donor: Inputs.getDonorDnaSample(meta),
            rna: Inputs.getTumorRnaSample(meta),
        ]

        meta_sample_entries.collect { sample_type, meta_sample ->
            if(meta_sample.containsKey(samplesheet.FileType.BAM){}

        }
    }



    // Sort inputs, separate by tumor and normal
    // channel: [ meta ]
    ch_inputs_tumor_sorted = ch_inputs
        .branch { meta ->
            def has_existing = Inputs.hasExistingInput(meta, Inputs.KEY.BAM_DNA_TUMOR)
            runnable: Inputs.hasTumorDnaFastq(meta) && !has_existing
            skip: true
        }

    ch_inputs_normal_sorted = ch_inputs
        .branch { meta ->
            def has_existing = Inputs.hasExistingInput(meta, Inputs.KEY.BAM_DNA_NORMAL)
            runnable: Inputs.hasNormalDnaFastq(meta) && !has_existing
            skip: true
        }

    ch_inputs_donor_sorted = ch_inputs
        .branch { meta ->
            def has_existing = Inputs.hasExistingInput(meta, Inputs.KEY.BAM_DNA_DONOR)
            runnable: Inputs.hasDonorDnaFastq(meta) && !has_existing
            skip: true
        }

    ch_inputs_rna_sorted = ch_inputs
        .branch { meta ->
            def has_existing = Inputs.hasExistingInput(meta, Inputs.KEY.BAM_RNA_TUMOR)
            runnable: Inputs.hasTumorRnaFastq(meta) && !has_existing
            skip: true
        }

    // Create FASTQ input channel
    // channel: [ meta_fastq, fastq_fwd, fastq_rev ]
    ch_fastq_inputs = Channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable.map { meta -> [meta, Inputs.getTumorDnaSample(meta), 'tumor'] },
            ch_inputs_normal_sorted.runnable.map { meta -> [meta, Inputs.getNormalDnaSample(meta), 'normal'] },
            ch_inputs_donor_sorted.runnable.map { meta -> [meta, Inputs.getDonorDnaSample(meta), 'donor'] },
            ch_inputs_rna_sorted.runnable.map { meta -> [meta, Inputs.getTumorRnaSample(meta), 'rna'] },
        )
        .flatMap { meta, meta_sample, sample_type ->
            meta_sample
                .getAt(samplesheet.FileType.FASTQ)
                .collect { key, fps ->
                    def (library_id, lane) = key

                    def sample_id = meta_sample.getOrDefault('longitudinal_sample_id', meta_sample['sample_id'])

                    def meta_fastq = [
                        key: meta.group_id,
                        id: "${meta.group_id}_${sample_id}",
                        sample_id: sample_id,
                        library_id: library_id,
                        lane: lane,
                        sample_type: sample_type,
                    ]

                    return [meta_fastq, fps['fwd'], fps['rev']]
                }
        }

    //
    // MODULE: fastp, fastq-tools
    //
    // UMI processing
    ch_fastq_umi_processed = Channel.empty()
    if(umi_enable) {

        def maybe_supported_panel = RefData.SupportedPanel.fromString(params.panel)

        if (maybe_supported_panel == RefData.SupportedPanel.MSK) {
            FASTQ_TOOLS_UMI_PROCESSING(
                ch_fastq_inputs,
                known_umis,
            )

            ch_fastq_umi_processed = FASTQ_TOOLS_UMI_PROCESSING.out.fastq
            ch_versions = ch_versions.mix(FASTQ_TOOLS_UMI_PROCESSING.out.versions)

        } else {
            FASTP_UMI_PROCESSING(
                ch_fastq_inputs,
                -1 // max_fastq_records
                umi_location,
                umi_length,
                umi_skip,
            )

            ch_fastq_umi_processed = FASTP_UMI_PROCESSING.out.fastq
            ch_versions = ch_versions.mix(FASTP_UMI_PROCESSING.out.versions)
        }

    } else {
        ch_fastq_umi_processed = ch_fastq_inputs
    }

    //
    // MODULE: fastp
    //
    // Split FASTQ into chunks if requested for distributed processing
    ch_fastqs_ready = Channel.empty()
    if (max_fastq_records > 0) {

        FASTP_FASTQ_SPLITTING(
            ch_fastq_umi_processed,
            max_fastq_records,
            "",
            0,
            -1,
        )

        ch_versions = ch_versions.mix(FASTP_UMI_PROCESSING.out.versions)

        ch_fastqs_ready = FASTP_FASTQ_SPLITTING.out.fastq
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
        ch_fastq_source = umi_enable ? FASTP_FASTQ_SPLITTING.out.fastq : ch_fastq_inputs

        ch_fastqs_ready = ch_fastq_source
            .map { meta_fastq, fastq_fwd, fastq_rev ->

                def meta_fastq_ready = [
                    *:meta_fastq,
                    split: null,
                ]

                return [meta_fastq_ready, fastq_fwd, fastq_rev]
            }

        ch_fastqs_ready = ch_fastq_umi_processed
    }

//     // Now prepare according to FASTQs splitting
//     // channel: [ meta_fastq_ready, fastq_fwd, fastq_fwd ]
//     ch_fastqs_ready = Channel.empty()
//     if (max_fastq_records > 0) {
//
//         ch_fastqs_ready = FASTP_FASTQ_SPLITTING.out.fastq
//             .flatMap { meta_fastq, reads_fwd, reads_rev ->
//
//                 def data = [reads_fwd, reads_rev]
//                     .transpose()
//                     .collect { fwd, rev ->
//
//                         def split_fwd = fwd.name.replaceAll('\\..+$', '')
//                         def split_rev = rev.name.replaceAll('\\..+$', '')
//
//                         assert split_fwd == split_rev
//
//                         // NOTE(SW): split allows meta_fastq_ready to be unique, which is required during reunite below
//                         def meta_fastq_ready = [
//                             *:meta_fastq,
//                             id: "${meta_fastq.id}_${split_fwd}",
//                             split: split_fwd,
//                         ]
//
//                         return [meta_fastq_ready, fwd, rev]
//                     }
//
//                 return data
//             }
//
//     } else {
//
//         // Select appropriate source
//         ch_fastq_source = umi_enable ? FASTP_FASTQ_SPLITTING.out.fastq : ch_fastq_inputs
//
//         ch_fastqs_ready = ch_fastq_source
//             .map { meta_fastq, fastq_fwd, fastq_rev ->
//
//                 def meta_fastq_ready = [
//                     *:meta_fastq,
//                     split: null,
//                 ]
//
//                 return [meta_fastq_ready, fastq_fwd, fastq_rev]
//             }
//
//     }

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
            channels.WorkflowChannels.restoreMeta(ch_bams_united.tumor, ch_inputs),
            channels.PlaceholderChannels.bamBai(ch_inputs_tumor_sorted.skip),
        )

    ch_bam_normal_out = Channel.empty()
        .mix(
            channels.WorkflowChannels.restoreMeta(ch_bams_united.normal, ch_inputs),
            channels.PlaceholderChannels.bamBai(ch_inputs_normal_sorted.skip),
        )

    ch_bam_donor_out = Channel.empty()
        .mix(
            channels.WorkflowChannels.restoreMeta(ch_bams_united.donor, ch_inputs),
            channels.PlaceholderChannels.bamBai(ch_inputs_donor_sorted.skip),
        )

    emit:
    dna_tumor  = ch_bam_tumor_out  // channel: [ meta, [bam, ...], [bai, ...] ]
    dna_normal = ch_bam_normal_out // channel: [ meta, [bam, ...], [bai, ...] ]
    dna_donor  = ch_bam_donor_out  // channel: [ meta, [bam, ...], [bai, ...] ]

    versions   = ch_versions       // channel: [ versions.yml ]
}
