//
// Apply post-alignment processing
//

import Constants
import Utils

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
    sequencing_platform    // string:  [mandatory] sequencing platform
    targeted_mode          // boolean: [mandatory] Set targeted mode
    umi_enable             // boolean: [mandatory] enable UMI processing
    umi_duplex_delim       // string:  [optional] UMI duplex delimiter
    generate_tsvs_only     // boolean: [mandatory] Generate REDUX TSVs from existing REDUX BAMs
    targeted_mode          // boolean: [mandatory] Set targeted mode

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources then sort, separating by sample type
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
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.REDUX_DIR_TUMOR)
            runnable: bams && ! has_existing
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
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.REDUX_DIR_NORMAL)
            runnable: bams && ! has_existing
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
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.REDUX_DIR_DONOR)
            runnable: bams && ! has_existing
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

            def sample_id = meta_sample.getOrDefault('longitudinal_sample_id', meta_sample['sample_id'])

            def meta_redux = [
                key: meta.group_id,
                id: "${meta.group_id}_${sample_id}",
                sample_id: sample_id,
                sample_type: sample_type,
            ]

            sample_data: [meta_redux, bams, bais]
            generate_tsvs_only: meta_sample.getOrDefault(Constants.InfoField.GENERATE_REDUX_TSVS_ONLY, false)
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
        sequencing_platform,
        targeted_mode,
        ch_redux_inputs.generate_tsvs_only,
        umi_enable,
        umi_duplex_delim,
        generate_tsvs_only,
        targeted_mode,
    )

    ch_versions = ch_versions.mix(REDUX.out.versions)

    // Split into sample type channels
    // channel: [ meta_redux, redux_dir ]
    ch_redux_out_sorted = REDUX.out.redux_dir
        .branch { meta_redux, redux_dir ->
            assert ['tumor', 'normal', 'donor'].contains(meta_redux.sample_type)
            tumor: meta_redux.sample_type == 'tumor'
            normal: meta_redux.sample_type == 'normal'
            donor: meta_redux.sample_type == 'donor'
            placeholder: true
        }

    // Set outputs, restoring original meta
    // channel: [ meta, redux_dir ]
    ch_outputs_tumor = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_redux_out_sorted.tumor, ch_inputs),
            ch_inputs_tumor_sorted.skip.map { meta -> [meta, []] },
        )

    // channel: [ meta, redux_dir ]
    ch_outputs_normal = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_redux_out_sorted.normal, ch_inputs),
            ch_inputs_normal_sorted.skip.map { meta -> [meta, []] },
        )

    // channel: [ meta, redux_dir ]
    ch_outputs_donor = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_redux_out_sorted.donor, ch_inputs),
            ch_inputs_donor_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    tumor_dir  = ch_outputs_tumor  // channel: [ meta, redux_dir ]
    normal_dir = ch_outputs_normal // channel: [ meta, redux_dir ]
    donor_dir  = ch_outputs_donor  // channel: [ meta, redux_dir ]

    versions   = ch_versions       // channel: [ versions.yml ]
}
