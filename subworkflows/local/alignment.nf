include { BWA_MEM        } from '../../modules/local/bwa/mem/main'
include { MARKDUPS       } from '../../modules/local/markdups/main'
include { FASTP          } from '../../modules/local/fastp/main'
include { SAMBAMBA_INDEX } from '../../modules/local/sambamba/index/main'

workflow ALIGNMENT {
    take:
    ch_inputs         // channel: [ meta ]
    genome_fasta
    genome_fai
    genome_dict
    genome_bwa_index
    unmap_regions
    max_fastq_lines

    main:
    // channel: [ group_id, sample_count ]
    ch_sample_counts = ch_inputs.map { meta -> [meta.group_id, Utils.groupSampleCounts(meta)] }

    // channel: [ meta ] (One sample per record).
    ch_meta_samples = ch_inputs.flatMap { meta -> Utils.splitGroupIntoSamples(meta) }

    // Sort inputs
    // channel: [ meta ] (One sample per record).
    ch_meta_samples_sorted = ch_meta_samples
        .branch { meta ->
            runnable_fastq: Utils.hasDnaFastq(meta)
            runnable_markdups: Utils.hasDnaMarkdupsBam(meta)
            skip: true
        }

    // channel: [ meta_fastq, reads_fwd_fastq, reads_rev_fastq ]
    ch_fastq_pairs = ch_meta_samples_sorted.runnable_fastq
        .flatMap { meta ->

            def sample_key = Constants.DNA_SAMPLE_KEYS.find { key -> meta.containsKey(key) }
            if (sample_key === null) {
                log.error "No DNA sample found"
                System.exit(1)
            }

            def sample_id = meta[sample_key]['sample_id']
            def fastq_files = meta[sample_key][Constants.FileType.FASTQ].tokenize(';')

            def meta_fastq_common = [:]
            meta.each { key, value ->


                if (key === sample_key) {
                    return
                }

                meta_fastq_common[key] = meta[key]
            }
            meta_fastq_common['sample_key'] = sample_key
            meta_fastq_common['sample_id'] = sample_id

            def fastq_pairs = []
            for (i = 0; i < fastq_files.size(); i += 2) {
                def reads_fwd = fastq_files[i]
                def reads_rev = fastq_files[i + 1]

                def meta_fastq = meta_fastq_common.getClass().newInstance(meta_fastq_common)
                meta_fastq['read_group'] = Utils.readGroupFromFastqPath(reads_fwd)

                fastq_pairs.add([meta_fastq, reads_fwd, reads_rev])
            }

            fastq_pairs
        }

    // Split fastq files using fastp.
    // channel: [ meta_fastq, reads_fwd_fastqs, reads_rev_fastqs ]
    ch_split_fastq_pairs = Channel.empty()
    if (max_fastq_lines > 0) {
        FASTP(ch_fastq_pairs)
        ch_split_fastq_pairs = FASTP.out.fastq
    } else {
        ch_split_fastq_pairs = ch_fastq_pairs.map { fastq_pair -> [fastq_pair[0], [fastq_pair[1]], [fastq_pair[2]]] }
    }

    // Create inputs for bwa mem.
    // channel: [ meta_fastq, reads_fwd_fastq, reads_rev_fastq ]
    ch_bwa_mem_inputs = ch_split_fastq_pairs.flatMap { fastq ->
        def meta = fastq[0]
        def fwd_reads = fastq[1]
        def rev_reads = fastq[2]

        // Pair up the reads.
        def read_pairs = [:]
        if (fwd_reads.size() == 1) {
            read_pairs[""] = ["000", fwd_reads[0], rev_reads[0]]
        } else {
            fwd_reads.each { fastq_path ->

                def base_name = fastq_path.getFileName().toString()
                def pattern = /^(\d+)\.(.+)_R[12]\.fastp\.fastq$/
                def matcher = base_name =~ pattern
                assert matcher.find()
                def split = matcher[0][1]
                def key = "${split}.${matcher[0][2]}"
                assert !read_pairs.containsKey(key)
                read_pairs[key] = [split, fastq_path]
            }

            rev_reads.each { fastq_path ->

                def base_name = fastq_path.getFileName().toString()
                def pattern = /^(.+)_R[12]\.fastp\.fastq$/
                def matcher = base_name =~ pattern
                assert matcher.find()
                def key = matcher[0][1]
                assert read_pairs.containsKey(key)
                read_pairs[key].add(fastq_path)
            }
        }

        def fastqs = []
        read_pairs.values().each { split_fastq_pair ->

            meta_fastq = meta.getClass().newInstance(meta)
            meta_fastq['split'] = split_fastq_pair[0]

            fastqs.add([meta_fastq, split_fastq_pair[1], split_fastq_pair[2]])
        }

        fastqs
    }

    // channel: [ meta_fastq, bam ]
    BWA_MEM(
        ch_bwa_mem_inputs,
        genome_fasta,
        genome_bwa_index,
    )

    // channel: [ meta_fastq, bam, bai ]
    SAMBAMBA_INDEX(
        BWA_MEM.out.bam,
    )

    // Prepare input to markdups process.
    // First we prepare a channel of inputs that have gone through alignment.
    // channel: [ meta_bam, bams, bais ]
    ch_fastq_markdups_inputs = SAMBAMBA_INDEX.out.bam
        .map { bam ->  // Strip read groups and splits.

            def meta = bam[0]
            def meta_bam = [:]
            meta.keySet().each { key ->

                if (key == 'read_group' || key == 'split') {
                    return
                }

                meta_bam[key] = meta[key]
            }

            [meta_bam, [meta_bam, bam[1], bam[2]]]
        }
        .groupTuple()
        .map { key_lane_bams ->
            def lane_bams = key_lane_bams[1]
            def meta_bam = lane_bams[0][0]
            def bams = []
            def bais = []
            lane_bams.each { lane_bam ->

                bams.add(lane_bam[1])
                bais.add(lane_bam[2])
            }

            [meta_bam, bams, bais]
        }

    // Next we prepare channel for markdups input that started of as aligned bams.
    // channel: [ meta, bams, bais ] (One sample per meta record).
    ch_input_markdups_inputs = ch_meta_samples_sorted.runnable_markdups.map { meta ->

        def sample_key = Constants.DNA_SAMPLE_KEYS.find { key -> meta.containsKey(key) }
        if (sample_key === null) {
            log.error "No DNA sample found"
            System.exit(1)
        }

        def sample_id = meta[sample_key]['sample_id']
        def bam = meta[sample_key][Constants.FileType.BAM_MARKDUPS]
        def bai = meta[sample_key][Constants.FileType.BAI_MARKDUPS]

        def meta_bam = meta.getClass().newInstance(meta);
        meta_bam['sample_key'] = sample_key
        meta_bam['sample_id'] = sample_id

        [meta_bam, [bam], [bai]]
    }

    // Merging the two markdups input channels.
    // channel: [ meta_bam, bams, bais ]
    ch_markdups_inputs = Channel.empty()
        .mix(
            ch_fastq_markdups_inputs,
            ch_input_markdups_inputs,
        )

    // channel: [ meta_bam, bam, bai ]
    MARKDUPS(
        ch_markdups_inputs,
        genome_fasta,
        genome_fai,
        genome_dict,
        unmap_regions,
    )

    // Fill the sample information back in.
    // channel: [ meta ] (One sample per meta record).
    ch_bam_samples = MARKDUPS.out.bam.map { bam ->

        def meta_bam = bam[0]

        def meta = meta_bam.getClass().newInstance(meta_bam)
        meta.remove('sample_key')
        meta.remove('sample_id')

        def sample = [sample_id: meta_bam.sample_id]
        sample[Constants.FileType.BAM] = bam[1]
        sample[Constants.FileType.BAI] = bam[2]
        meta[meta_bam.sample_key] = sample

        meta
    }

    // Merge back in skipped meta entries.
    // channel: [ meta ] (One sample per meta record).
    ch_all_samples = Channel.empty()
        .mix(
            ch_bam_samples,
            ch_meta_samples_sorted.skip,
        )

    // Merge individual sample records back into group records without blocking for the whole channel to be processed.
    // channel: [ meta_bam ]
    ch_outputs = ch_sample_counts
        .cross(
            ch_all_samples.map { meta -> [meta.group_id, meta] }
        )
        .map { count_tuple, meta_tuple ->
            def group_id = count_tuple[0]
            def count = count_tuple[1]
            def meta = meta_tuple[1]

            tuple(groupKey(group_id, count), meta)
        }
        .groupTuple()
        .map { key_samples ->

            def samples = key_samples[1]
            def merged_sample = [:]
            samples.each { sample ->

                sample.each { key, value -> merged_sample[key] = value }
            }

            merged_sample
        }

    emit:
    meta_bam = ch_outputs
    // TODO[MC]: Channel version outputs.
}
