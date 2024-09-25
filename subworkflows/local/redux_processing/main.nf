//
// Apply post-alignment processing
//

import Constants
import Utils

include { REDUX } from '../../../modules/local/redux/main'

workflow REDUX_PROCESSING {
    take:
    // Sample data
    ch_inputs     // channel: [mandatory] [ meta ]
    ch_dna_tumor  // channel: [mandatory] [ meta, [bam, ...], [bai, ...] ]
    ch_dna_normal // channel: [mandatory] [ meta, [bam, ...], [bai, ...] ]
    ch_dna_donor  // channel: [optional]  [ meta, [bam, ...], [bai, ...] ]

    // Reference data
    genome_fasta     // channel: [mandatory] /path/to/genome_fasta
    genome_ver       // channel: [mandatory] genome version
    genome_fai       // channel: [mandatory] /path/to/genome_fai
    genome_dict      // channel: [mandatory] /path/to/genome_dict
    unmap_regions    // channel: [mandatory] /path/to/unmap_regions
    msi_jitter_sites // channel: [mandatory] /path/to/msi_jitter_sites

    // Params
    has_umis         // boolean: [mandatory] UMI processing flag
    umi_duplex_delim // string:  [optional] UMI duplex delimiter

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select and sort input sources, separating bytumor and normal
    // channel: runnable: [ meta, [bam, ...], [bai, ...] ]
    // channel: skip: [ meta ]
    ch_inputs_tumor_sorted = ch_dna_tumor
        .map { meta, bams, bais ->
            return [
                meta,
                Utils.hasExistingInput(meta, Constants.INPUT.BAM_DNA_TUMOR) ? [Utils.getInput(meta, Constants.INPUT.BAM_DNA_TUMOR)] : bams,
                Utils.hasExistingInput(meta, Constants.INPUT.BAI_DNA_TUMOR) ? [Utils.getInput(meta, Constants.INPUT.BAI_DNA_TUMOR)] : bais,
            ]
        }
        .branch { meta, bams, bais ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAM_REDUX_DNA_TUMOR)
            runnable: bams && !has_existing
            skip: true
                return meta
        }

    ch_inputs_normal_sorted = ch_dna_normal
        .map { meta, bams, bais ->
            return [
                meta,
                Utils.hasExistingInput(meta, Constants.INPUT.BAM_DNA_NORMAL) ? [Utils.getInput(meta, Constants.INPUT.BAM_DNA_NORMAL)] : bams,
                Utils.hasExistingInput(meta, Constants.INPUT.BAI_DNA_NORMAL) ? [Utils.getInput(meta, Constants.INPUT.BAI_DNA_NORMAL)] : bais,
            ]
        }
        .branch { meta, bams, bais ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAM_REDUX_DNA_NORMAL)
            runnable: bams && !has_existing
            skip: true
                return meta
        }

    ch_inputs_donor_sorted = ch_dna_donor
        .map { meta, bams, bais ->
            return [
                meta,
                Utils.hasExistingInput(meta, Constants.INPUT.BAM_DNA_DONOR) ? [Utils.getInput(meta, Constants.INPUT.BAM_DNA_DONOR)] : bams,
                Utils.hasExistingInput(meta, Constants.INPUT.BAI_DNA_DONOR) ? [Utils.getInput(meta, Constants.INPUT.BAI_DNA_DONOR)] : bais,
            ]
        }
        .branch { meta, bams, bais ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAM_REDUX_DNA_DONOR)
            runnable: bams && !has_existing
            skip: true
            return meta
        }

    // Create process input channel
    // channel: [ meta_redux, [bam, ...], [bai, ...] ]
    ch_redux_inputs = Channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable.map { meta, bams, bais -> [meta, Utils.getTumorDnaSample(meta), 'tumor', bams, bais] },
            ch_inputs_normal_sorted.runnable.map { meta, bams, bais -> [meta, Utils.getNormalDnaSample(meta), 'normal', bams, bais] },
            ch_inputs_donor_sorted.runnable.map { meta, bams, bais -> [meta, Utils.getDonorDnaSample(meta), 'donor', bams, bais] },
        )
        .map { meta, meta_sample, sample_type, bams, bais ->

            def meta_redux = [
                key: meta.group_id,
                id: "${meta.group_id}_${meta_sample.sample_id}",
                sample_id: meta_sample.sample_id,
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
        has_umis,
        umi_duplex_delim,
    )

    // Sort into a tumor and normal channel
    ch_redux_out = REDUX.out.bam
        .branch { meta, bam, bai ->
            assert ['tumor', 'normal', 'donor'].contains(meta.sample_type)
            tumor: meta.sample_type == 'tumor'
            normal: meta.sample_type == 'normal'
            donor: meta.sample_type == 'donor'
            placeholder: true
        }

    // Set outputs, restoring original meta
    // channel: [ meta, bam, bai ]
    ch_redux_tumor_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_redux_out.tumor, ch_inputs),
            ch_inputs_tumor_sorted.skip.map { meta -> [meta, [], []] },
        )

    ch_redux_normal_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_redux_out.normal, ch_inputs),
            ch_inputs_normal_sorted.skip.map { meta -> [meta, [], []] },
        )

    ch_redux_donor_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_redux_out.donor, ch_inputs),
            ch_inputs_donor_sorted.skip.map { meta -> [meta, [], []] },
        )

    emit:
    dna_tumor  = ch_redux_tumor_out  // channel: [ meta, bam, bai ]
    dna_normal = ch_redux_normal_out // channel: [ meta, bam, bai ]
    dna_donor = ch_redux_donor_out   // channel: [ meta, bam, bai ]

    versions   = ch_versions       // channel: [ versions.yml ]
}