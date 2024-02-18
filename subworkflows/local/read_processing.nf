include { MARKDUPS       } from '../../modules/local/markdups/main'

workflow READ_PROCESSING {
    take:
    // Sample data
    ch_inputs   // channel: [mandatory] [ meta ]
    ch_dna_bams // channel: [mandatory] [ meta, bam_dna ]
    ch_rna_bams // channel: [mandatory] [ meta, bam_rna ]
    genome_fasta
    genome_fai
    genome_dict
    unmap_regions

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // channel: [ group_id, sample_count ]
    ch_sample_counts = ch_inputs.map { meta -> [meta.group_id, Utils.groupSampleCounts(meta)] }

    // channel: [ meta ] (One sample per record).
    ch_meta_samples = ch_dna_bams.flatMap { meta -> Utils.splitGroupIntoSamples(meta) }

    // Sort inputs
    // channel: [ meta ] (One sample per record).
    ch_meta_samples_sorted = ch_meta_samples
        .branch { meta ->
            runnable: Utils.hasDnaMarkdupsBam(meta)
            skip: true
        }

    // MarkDups
    // Prepare input to markdups process.
    // channel: [ meta_bam, bams, bais ]
    ch_markdups_inputs = ch_meta_samples_sorted.runnable
        .map { meta_sample ->

            def meta_bam = Utils.shallow_copy(meta_sample)
            def bams = []
            def bais = []
            meta_sample.each { key, value ->

                if ((value instanceof java.util.Map) && value.containsKey('sample_id')) {
                    meta_bam['sample_id'] = value.sample_id
                    meta_bam['sample_key'] = key
                    bams = value[Constants.FileType.BAM_MARKDUPS]
                    bais = value[Constants.FileType.BAI_MARKDUPS]
                }
            }

            if (!(bams instanceof Collection)) {
                bams = [bams]
            }

            if (!(bais instanceof Collection)) {
                bais = [bais]
            }

            [meta_bam, bams, bais]
        }

    // channel: [ meta_bam, bam, bai ]
    MARKDUPS(
        ch_markdups_inputs,
        genome_fasta,
        genome_fai,
        genome_dict,
        unmap_regions,
        false,
    )

    ch_versions = ch_versions.mix(MARKDUPS.out.versions)

    // Update sample information.
    // channel: [ meta ] (One sample per meta record).
    ch_bam_samples = MARKDUPS.out.bam.map { bam ->

        def meta_bam = bam[0]

        def meta = Utils.shallow_copy(meta_bam)
        meta.remove('sample_id')
        meta.remove('sample_key')

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
    ch_markduplicates_dna_out = ch_sample_counts
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

    // TODO(SW): implement outputs
    ch_markduplicates_rna_out = Channel.empty()

    emit:
    dna       = ch_markduplicates_dna_out // channel: [ meta ]
    rna       = ch_markduplicates_rna_out // channel: [ meta, bam_rna ]
    versions  = ch_versions               // channel: [ versions.yml ]
}
