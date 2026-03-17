//
// Apply post-alignment processing
//

include { REDUX } from '../../../modules/local/redux/main'

workflow REDUX_PROCESSING {
    take:
    // Sample data
    ch_inputs        // channel: [mandatory] [ meta ]
    ch_dna_tumor     // channel: [mandatory] [ meta, [bam, ...], [bai, ...] ]
    ch_dna_normal    // channel: [mandatory] [ meta, [bam, ...], [bai, ...] ]
    ch_dna_donor     // channel: [mandatory] [ meta, [bam, ...], [bai, ...] ]

    // Reference data
    genome_fasta     // channel: [mandatory] /path/to/genome_fasta
    genome_ver       // channel: [mandatory] genome version
    genome_fai       // channel: [mandatory] /path/to/genome_fai
    genome_dict      // channel: [mandatory] /path/to/genome_dict
    unmap_regions    // channel: [mandatory] /path/to/unmap_regions
    msi_jitter_sites // channel: [mandatory] /path/to/msi_jitter_sites

    // Params
    sequencing_type  // string:  [mandatory] sequencing type
    umi_enable       // boolean: [mandatory] enable UMI processing
    umi_duplex_delim // string:  [optional] UMI duplex delimiter
    targeted_mode    // boolean: [mandatory] Set targeted mode

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select and sort input sources, separating by sample type
    // channel: runnable: [ meta, [bam, ...], [bai, ...] ]
    // channel: skip: [ meta ]

    def selectBamInputs = { ch_dna, bam_type, bai_type, bam_redux_type ->
        return ch_dna.map { meta, bams, bais ->

            bams = Inputs.hasExistingInput(meta, bam_type)
                ? [Inputs.getInput(meta, bam_type)]
                : bams

            bais = Inputs.hasExistingInput(meta, bai_type)
                ? [Inputs.getInput(meta, bai_type)]
                : bais

            return [meta, bams, bais]
        }
        .branch { meta, bams, bais ->
            def has_existing = Inputs.hasExistingInput(meta, bam_redux_type)
            runnable: bams && !has_existing
            skip: true
                return meta
        }
    }

    ch_inputs_tumor = selectBamInputs(
        ch_dna_tumor,
        SampleMeta.INPUT.BAM_DNA_TUMOR,
        SampleMeta.INPUT.BAI_DNA_TUMOR,
        SampleMeta.INPUT.BAM_REDUX_DNA_TUMOR
    )

    ch_inputs_normal = selectBamInputs(
        ch_dna_normal,
        SampleMeta.INPUT.BAM_DNA_NORMAL,
        SampleMeta.INPUT.BAI_DNA_NORMAL,
        SampleMeta.INPUT.BAM_REDUX_DNA_NORMAL
    )

    ch_inputs_donor = selectBamInputs(
        ch_dna_donor,
        SampleMeta.INPUT.BAM_DNA_DONOR,
        SampleMeta.INPUT.BAI_DNA_DONOR,
        SampleMeta.INPUT.BAM_REDUX_DNA_DONOR
    )

    // Create process input channel
    // channel: [ meta_redux, [bam, ...], [bai, ...] ]
    ch_redux_inputs = Channel.empty()
        .mix(
            ch_inputs_tumor.runnable.map { meta, bams, bais -> [meta, Inputs.getTumorDnaSample(meta), 'tumor', bams, bais] },
            ch_inputs_normal.runnable.map { meta, bams, bais -> [meta, Inputs.getNormalDnaSample(meta), 'normal', bams, bais] },
            ch_inputs_donor.runnable.map { meta, bams, bais -> [meta, Inputs.getDonorDnaSample(meta), 'donor', bams, bais] },
        )
        .map { meta, meta_sample, sample_type, bams, bais ->

            def sample_id = meta_sample.getOrDefault('longitudinal_sample_id', meta_sample['sample_id'])

            def meta_redux = [
                key: meta.group_id,
                id: "${meta.group_id}_${sample_id}",
                sample_id: sample_id,
                sample_type: sample_type,
            ]

            return [meta_redux, bams, bais]
        }

    // Run process
    REDUX(
        ch_redux_inputs,
        genome_fasta,
        genome_ver,
        genome_fai,
        genome_dict,
        unmap_regions,
        msi_jitter_sites,
        sequencing_type,
        umi_enable,
        umi_duplex_delim,
        targeted_mode,
    )

    ch_versions = ch_versions.mix(REDUX.out.versions)

    // Combine TSV outputs into single channel for processing
    // channel: [ meta, bam, bai, bqr_tsv, jitter_tsv, ms_tsv, bqr_plot ]
    ch_redux_out = WorkflowOncoanalyser.groupByMeta(
        REDUX.out.bam,
        REDUX.out.bqr_tsv,
        REDUX.out.dup_freq_tsv,
        REDUX.out.jitter_tsv,
        REDUX.out.ms_tsv,
        REDUX.out.bqr_plot,
    )

    // Sort into a tumor and normal channel
    // channel: [ meta, bam, bai, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv, bqr_plot ]
    ch_redux_out_sorted = ch_redux_out
        .branch { meta, bam, bai, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv, bqr_plot ->
            assert ['tumor', 'normal', 'donor'].contains(meta.sample_type)
            tumor: meta.sample_type == 'tumor'
            normal: meta.sample_type == 'normal'
            donor: meta.sample_type == 'donor'
            placeholder: true
        }

    // Set outputs, restoring original meta, split by file type
    // channel: [ meta, bam, bai, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv, bqr_plot ]
    def createOutputChannels = { ch_redux_out_sample_type, ch_sample_type_skip ->

        def placeholder_bam = [[]] * PlaceholderChannels.N_ITEMS_BAM_BAI
        def placeholder_tsv = [[]] * PlaceholderChannels.N_ITEMS_REDUX_TSVS
        def placeholder_plot = [[]] * PlaceholderChannels.N_ITEMS_REDUX_PLOTS
        def placeholders = [*placeholder_bam, *placeholder_tsv, *placeholder_plot]

        return Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(ch_redux_out_sample_type, ch_inputs),
                ch_sample_type_skip.map { meta -> [meta, *placeholders] },
            )
            .multiMap { meta, bam, bai, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv, bqr_plot ->
                bam: [meta, bam, bai]
                tsv: [meta, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv]
                plot: [meta, bqr_plot]
            }
    }

    ch_redux_tumor_out = createOutputChannels(ch_redux_out_sorted.tumor, ch_inputs_tumor.skip)
    ch_redux_normal_out = createOutputChannels(ch_redux_out_sorted.normal, ch_inputs_normal.skip)
    ch_redux_donor_out = createOutputChannels(ch_redux_out_sorted.donor, ch_inputs_donor.skip)

    emit:
    dna_tumor_bam  = ch_redux_tumor_out.bam  // channel: [ meta, bam, bai ]
    dna_normal_bam = ch_redux_normal_out.bam // channel: [ meta, bam, bai ]
    dna_donor_bam  = ch_redux_donor_out.bam  // channel: [ meta, bam, bai ]

    dna_tumor_tsv  = ch_redux_tumor_out.tsv  // channel: [ meta, redux_tsv, ... ]
    dna_normal_tsv = ch_redux_normal_out.tsv // channel: [ meta, redux_tsv, ... ]
    dna_donor_tsv  = ch_redux_donor_out.tsv  // channel: [ meta, redux_tsv, ... ]

    dna_tumor_plot  = ch_redux_tumor_out.plot  // channel: [ meta, bqr_plot ]
    dna_normal_plot = ch_redux_normal_out.plot // channel: [ meta, bqr_plot ]
    dna_donor_plot  = ch_redux_donor_out.plot  // channel: [ meta, bqr_plot ]

    versions       = ch_versions             // channel: [ versions.yml ]
}
