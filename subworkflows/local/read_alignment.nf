include { BWA_MEM2 } from '../../modules/local/bwa/mem2/main'
include { FASTP          } from '../../modules/local/fastp/main'
include { SAMBAMBA_INDEX } from '../../modules/local/sambamba/index/main'
include { STAR     } from '../../modules/local/star/main'

workflow READ_ALIGNMENT {
    take:
    // Sample data
    ch_inputs // channel: [mandatory] [ meta ]
    genome_fasta
    genome_bwa_index
    max_fastq_records

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // channel: [ group_id, sample_count ]
    ch_sample_counts = ch_inputs.map { meta -> [meta.group_id, Utils.groupSampleCounts(meta)] }

    // channel: [ meta ] (One sample per record).
    ch_meta_samples = ch_inputs.flatMap { meta -> Utils.splitGroupIntoSamples(meta) }

    // Sort inputs
    // channel: [ meta ] (One sample per record).
    ch_meta_samples_sorted = ch_meta_samples
        .branch { meta ->
            runnable_fastq: Utils.hasDnaFastq(meta)
            skip: true
        }

    // STAR
    // TODO(SW): implement inputs
    // ch_star_inputs = Channel.of([[id: 'foo'], []])
    // STAR(
    //     ch_star_inputs,
    //     // TODO(SW): include reference files
    // )
    // TODO(SW): implement outputs
    ch_star_outputs = Channel.empty()

    // BWA MEM2
    // channel: [ sample_key, fastq_pair_count ]
    ch_sample_fastq_pair_count = ch_meta_samples_sorted.runnable_fastq.map { meta_sample ->

        def sample_key = Utils.shallow_copy(meta_sample)
        def fastq_pair_count = 0
        meta_sample.each { key, value ->

            if ((value instanceof java.util.Map) && value.containsKey('sample_id')) {
                sample_key['sample_key'] = key
                sample_key['sample_id'] = value.sample_id
                sample_key.remove(key)

                fastq_pair_count = value[Constants.FileType.FASTQ].toString().tokenize(';').size() / 2
            }
        }

        [sample_key, fastq_pair_count]
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
            def fastq_files = meta[sample_key][Constants.FileType.FASTQ].toString().tokenize(';')

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

                def meta_fastq = Utils.shallow_copy(meta_fastq_common)
                meta_fastq['read_group'] = Utils.readGroupFromFastqPath(reads_fwd)

                fastq_pairs.add([meta_fastq, reads_fwd, reads_rev])
            }

            fastq_pairs
        }

    // Split fastq files using fastp.
    // channel: [ meta_fastq, reads_fwd_fastqs, reads_rev_fastqs ]
    ch_split_fastq_pairs = Channel.empty()
    if (max_fastq_records > 0) {
        FASTP(
            ch_fastq_pairs,
            max_fastq_records,
        )

        ch_versions = ch_versions.mix(FASTP.out.versions)

        ch_split_fastq_pairs = FASTP.out.fastq
    } else {
        ch_split_fastq_pairs = ch_fastq_pairs.map { fastq_pair -> [fastq_pair[0], [fastq_pair[1]], [fastq_pair[2]]] }
    }

    // channel: [ sample_key, fastq_pair_split_count ]
    ch_sample_fastq_pair_split_count = ch_sample_fastq_pair_count
        .cross(
            ch_split_fastq_pairs.map { split_fastq_pairs ->

                def meta_sample = split_fastq_pairs[0]
                def sample_key = Utils.shallow_copy(meta_sample)
                sample_key.remove('read_group')
                sample_key.remove(meta_sample.sample_key)

                [sample_key, split_fastq_pairs[1].size()]
            }
        )
        .map { count_tuple, split_count_tuple ->
            def sample_key = count_tuple[0]
            def count = count_tuple[1].intValue()
            def split_count = split_count_tuple[1]

            tuple(groupKey(sample_key, count), sample_key, split_count)
        }
        .groupTuple()
        .map { group_key, sample_keys, split_counts ->

            [sample_keys[0], split_counts.sum()]
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

            meta_fastq = Utils.shallow_copy(meta)
            meta_fastq['split'] = split_fastq_pair[0]

            fastqs.add([meta_fastq, split_fastq_pair[1], split_fastq_pair[2]])
        }

        fastqs
    }

    // channel: [ meta_fastq, bam ]
    BWA_MEM2(
        ch_bwa_mem_inputs,
        genome_fasta,
        genome_bwa_index,
    )

    ch_versions = ch_versions.mix(BWA_MEM2.out.versions)

    // channel: [ meta_fastq, bam, bai ]
    SAMBAMBA_INDEX(
        BWA_MEM2.out.bam,
    )

    ch_versions = ch_versions.mix(SAMBAMBA_INDEX.out.versions)

    // Merge all bam records for a single sample into a singlke record.
    // channel: [ meta ] (One sample per meta record).
    ch_merged_bam_samples = ch_sample_fastq_pair_split_count
        .cross(
            SAMBAMBA_INDEX.out.bam
                .map { bam ->

                    def meta_bam = bam[0]
                    def sample_key = Utils.shallow_copy(meta_bam)
                    sample_key.remove(meta_bam.sample_key)
                    sample_key.remove('read_group')
                    sample_key.remove('split')

                    [sample_key, bam]
                }
        )
        .map { count_tuple, bam_tuple ->

            def sample_key = count_tuple[0]
            def count = count_tuple[1]
            def bam = bam_tuple[1]

            tuple(groupKey(sample_key, count), bam)
        }
        .groupTuple()
        .map { group_key, bams ->

            def first_meta_bam = bams[0][0]
            def sample_key = first_meta_bam.sample_key

            def bam_files = []
            def bai_files = []

            def meta_bam = Utils.shallow_copy(first_meta_bam)
            meta_bam.remove(sample_key)
            meta_bam.remove('sample_key')
            meta_bam.remove('sample_id')
            meta_bam.remove('read_group')
            meta_bam.remove('split')

            meta_bam[sample_key] = [sample_id: first_meta_bam.sample_id]
            meta_bam[sample_key][Constants.FileType.BAM_MARKDUPS] = bam_files
            meta_bam[sample_key][Constants.FileType.BAI_MARKDUPS] = bai_files

            bams.each { bam ->

                bam_files.add(bam[1])
                bai_files.add(bam[2])
            }

            meta_bam
        }

    // Merge back in skipped meta entries.
    // channel: [ meta ] (One sample per meta record).
    ch_all_samples = Channel.empty()
        .mix(
            ch_merged_bam_samples,
            ch_meta_samples_sorted.skip,
        )

    // Merge individual sample records back into group records without blocking for the whole channel to be processed.
    // channel: [ meta_bam ]
    ch_bwa_outputs = ch_sample_counts
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
        .map { group_key, meta_samples ->

            def meta_group = [:]
            meta_samples.each { meta_sample ->

                meta_sample.each { key, value -> meta_group[key] = value }
            }

            meta_group
        }

    emit:
    dna       = ch_bwa_outputs  // channel: [ meta ]
    // TODO(SW): RNA alignment.
    rna       = ch_star_outputs // channel: [ meta, bam_rna ]
    versions  = ch_versions     // channel: [ versions.yml ]
}
