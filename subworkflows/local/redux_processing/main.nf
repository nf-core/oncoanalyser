//
// Apply post-alignment processing
//

import Constants
import Utils

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

            bams = Utils.hasExistingInput(meta, bam_type)
                ? [Utils.getInput(meta, bam_type)]
                : bams

            bais = Utils.hasExistingInput(meta, bai_type)
                ? [Utils.getInput(meta, bai_type)]
                : bais

            return [meta, bams, bais]
        }
        .branch { meta, bams, bais ->
            def has_existing = Utils.hasExistingInput(meta, bam_redux_type)
            runnable: bams && !has_existing
            skip: true
                return meta
        }
    }

    ch_inputs_tumor = selectBamInputs(
        ch_dna_tumor,
        Constants.INPUT.BAM_DNA_TUMOR,
        Constants.INPUT.BAI_DNA_TUMOR,
        Constants.INPUT.BAM_REDUX_DNA_TUMOR
    )

    ch_inputs_normal = selectBamInputs(
        ch_dna_normal,
        Constants.INPUT.BAM_DNA_NORMAL,
        Constants.INPUT.BAI_DNA_NORMAL,
        Constants.INPUT.BAM_REDUX_DNA_NORMAL
    )

    ch_inputs_donor = selectBamInputs(
        ch_dna_donor,
        Constants.INPUT.BAM_DNA_DONOR,
        Constants.INPUT.BAI_DNA_DONOR,
        Constants.INPUT.BAM_REDUX_DNA_DONOR
    )

    // Create process input channel
    // channel: [ meta_redux, [bam, ...], [bai, ...] ]
    ch_redux_inputs = Channel.empty()
        .mix(
            ch_inputs_tumor.runnable.map { meta, bams, bais -> [meta, Utils.getTumorDnaSample(meta), 'tumor', bams, bais] },
            ch_inputs_normal.runnable.map { meta, bams, bais -> [meta, Utils.getNormalDnaSample(meta), 'normal', bams, bais] },
            ch_inputs_donor.runnable.map { meta, bams, bais -> [meta, Utils.getDonorDnaSample(meta), 'donor', bams, bais] },
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

    // Handle sequencing type
    // NOTE(LN): Comparing enums here at the workflow level allows the REDUX process to resume
    // by avoiding enum imports in the REDUX process that would lead to task serialisation issues
    sequencing_type_enum = Utils.getEnumFromString(sequencing_type, Constants.SequencingType)
    sequencing_type_ultima = sequencing_type_enum === Constants.SequencingType.ULTIMA

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
        sequencing_type_ultima,
        umi_enable,
        umi_duplex_delim,
        targeted_mode,
    )

    ch_versions = ch_versions.mix(REDUX.out.versions)

    // Combine TSV outputs into single channel for processing
    // channel: [ meta, [bam, bai], [bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv] ]
    ch_redux_out = WorkflowOncoanalyser.groupByMeta(
        ["flatten": false],
        REDUX.out.bam,
        REDUX.out.tsv,
    )

    // Sort into a tumor and normal channel
    // channel: [ meta, [bam, bai], [bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv] ]
    ch_redux_out_sorted = ch_redux_out
        .branch { meta, bam, tsv ->
            assert ['tumor', 'normal', 'donor'].contains(meta.sample_type)
            tumor: meta.sample_type == 'tumor'
            normal: meta.sample_type == 'normal'
            donor: meta.sample_type == 'donor'
            placeholder: true
        }

    // Set outputs, restoring original meta, split into BAMs and TSVs
    // channel: [ meta, bam, bai, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv ]
    def createOutputChannels = { ch_sample_type_redux_out, ch_sample_type_skip ->

        def empty_bam = [[],[]]
        def empty_tsv = [[],[],[],[]]

        return Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(ch_sample_type_redux_out, ch_inputs),
                ch_sample_type_skip.map { meta -> [meta, empty_bam, empty_tsv] },
            )
            .multiMap { meta, bam, tsv ->
                bam: [meta, *bam]
                tsv: [meta, *tsv]
            }
    }

    ch_redux_tumor_out = createOutputChannels(ch_redux_out_sorted.tumor, ch_inputs_tumor.skip)
    ch_redux_normal_out = createOutputChannels(ch_redux_out_sorted.normal, ch_inputs_normal.skip)
    ch_redux_donor_out = createOutputChannels(ch_redux_out_sorted.donor, ch_inputs_donor.skip)

    emit:
    dna_tumor_bam  = ch_redux_tumor_out.bam  // channel: [ meta, bam, bai ]
    dna_normal_bam = ch_redux_normal_out.bam // channel: [ meta, bam, bai ]
    dna_donor_bam  = ch_redux_donor_out.bam  // channel: [ meta, bam, bai ]

    dna_tumor_tsv  = ch_redux_tumor_out.tsv  // channel: [ meta, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv ]
    dna_normal_tsv = ch_redux_normal_out.tsv // channel: [ meta, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv ]
    dna_donor_tsv  = ch_redux_donor_out.tsv  // channel: [ meta, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv ]

    versions       = ch_versions             // channel: [ versions.yml ]
}
