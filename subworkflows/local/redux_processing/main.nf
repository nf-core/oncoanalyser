//
// Apply post-alignment processing
//

include { REDUX } from '../../../modules/local/redux/main'

workflow REDUX_PROCESSING {
    take:
    // Sample data
    ch_inputs              // channel: [mandatory] [ meta ]
    ch_dna_tumor           // channel: [mandatory] [ meta, [bam, ...], [bai, ...] ]
    ch_dna_normal          // channel: [mandatory] [ meta, [bam, ...], [bai, ...] ]
    ch_dna_donor           // channel: [mandatory] [ meta, [bam, ...], [bai, ...] ]

    // Reference data
    genome_fasta           // channel: [mandatory] /path/to/genome_fasta
    genome_ver             // channel: [mandatory] genome version
    genome_fai             // channel: [mandatory] /path/to/genome_fai
    genome_dict            // channel: [mandatory] /path/to/genome_dict
    unmap_regions          // channel: [mandatory] /path/to/unmap_regions
    msi_jitter_sites       // channel: [mandatory] /path/to/msi_jitter_sites
    msi_model_coefficients // channel: [mandatory] /path/to/msi_model_coefficients
    msi_model_error_rates  // channel: [mandatory] /path/to/msi_model_error_rates

    // Params
    sequencing_type        // string:  [mandatory] sequencing type
    umi_enable             // boolean: [mandatory] enable UMI processing
    umi_duplex_delim       // string:  [optional] UMI duplex delimiter
    targeted_mode          // boolean: [mandatory] Set targeted mode

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
        Inputs.KEY.BAM_DNA_TUMOR,
        Inputs.KEY.BAI_DNA_TUMOR,
        Inputs.KEY.BAM_REDUX_DNA_TUMOR
    )

    ch_inputs_normal = selectBamInputs(
        ch_dna_normal,
        Inputs.KEY.BAM_DNA_NORMAL,
        Inputs.KEY.BAI_DNA_NORMAL,
        Inputs.KEY.BAM_REDUX_DNA_NORMAL
    )

    ch_inputs_donor = selectBamInputs(
        ch_dna_donor,
        Inputs.KEY.BAM_DNA_DONOR,
        Inputs.KEY.BAI_DNA_DONOR,
        Inputs.KEY.BAM_REDUX_DNA_DONOR
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
        msi_model_coefficients,
        msi_model_error_rates,
        sequencing_type,
        umi_enable,
        umi_duplex_delim,
        targeted_mode,
    )

    ch_versions = ch_versions.mix(REDUX.out.versions)

    // Combine outputs into single channel for processing
    // channel: [ meta, bam, bai, dir ]
    ch_redux_out = channels.WorkflowChannels.groupByMeta(
        REDUX.out.bam,
        REDUX.out.dir,
    )

    // Split into sample type channels
    // channel: [ meta, bam, bai, dir ]
    ch_redux_out_sorted = ch_redux_out
        .branch { meta, bam, bai, dir ->
            assert ['tumor', 'normal', 'donor'].contains(meta.sample_type)
            tumor: meta.sample_type == 'tumor'
            normal: meta.sample_type == 'normal'
            donor: meta.sample_type == 'donor'
            placeholder: true
        }

    // Restore original meta, add placeholder channels for skipped samples, split by file type
    // channel: bam: [ meta, bam, bai ]
    // channel: dir: [ meta, dir ]
    def createOutputChannels = { ch_redux_out_sample_type, ch_sample_type_skip ->

        def placeholder_bam_bai = [[]] * channels.PlaceholderChannels.N_ITEMS_BAM_BAI
        def placeholder_dir = [[]] * channels.PlaceholderChannels.N_ITEMS_TOOL_DIR
        def ch_output_skip = ch_sample_type_skip.map { meta -> [meta, *placeholder_bam_bai, *placeholder_dir] }

        return Channel.empty()
            .mix(
                channels.WorkflowChannels.restoreMeta(ch_redux_out_sample_type, ch_inputs),
                ch_output_skip,
            )
            .multiMap { meta, bam, bai, dir ->
                bam: [meta, bam, bai]
                dir: [meta, dir]
            }
    }

    ch_redux_tumor_out = createOutputChannels(ch_redux_out_sorted.tumor, ch_inputs_tumor.skip)
    ch_redux_normal_out = createOutputChannels(ch_redux_out_sorted.normal, ch_inputs_normal.skip)
    ch_redux_donor_out = createOutputChannels(ch_redux_out_sorted.donor, ch_inputs_donor.skip)

    emit:
    dna_tumor_bam  = ch_redux_tumor_out.bam  // channel: [ meta, bam, bai ]
    dna_normal_bam = ch_redux_normal_out.bam // channel: [ meta, bam, bai ]
    dna_donor_bam  = ch_redux_donor_out.bam  // channel: [ meta, bam, bai ]

    dna_tumor_dir  = ch_redux_tumor_out.dir  // channel: [ meta, dir ]
    dna_normal_dir = ch_redux_normal_out.dir // channel: [ meta, dir ]
    dna_donor_dir  = ch_redux_donor_out.dir  // channel: [ meta, dir ]

    versions       = ch_versions             // channel: [ versions.yml ]
}
