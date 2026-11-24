//
// Apply post-alignment processing
//

include { REDUX } from '../../../modules/local/redux/main'

workflow REDUX_PROCESSING {
    take:
    // Sample data
    ch_inputs              // channel: [mandatory] [ meta ]
    ch_dna_tumor           // channel: [mandatory] [ meta, [aln, ...], [idx, ...] ]
    ch_dna_normal          // channel: [mandatory] [ meta, [aln, ...], [idx, ...] ]
    ch_dna_donor           // channel: [mandatory] [ meta, [aln, ...], [idx, ...] ]

    // Reference data
    genome_fasta           // channel: [mandatory] /path/to/genome_fasta
    genome_version         // channel: [mandatory] genome version
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

    main:
    // Select input sources then sort, separating by sample type
    // channel: runnable: [ meta, [aln, ...], [idx, ...] ]
    // channel: skip: [ meta ]
    ch_inputs_tumor_sorted = ch_dna_tumor
        .map { meta, alns, idxs ->
            return [
                meta,
                Utils.hasExistingInput(meta, Constants.INPUT.ALN_DNA_TUMOR) ? [Utils.getInput(meta, Constants.INPUT.ALN_DNA_TUMOR)] : alns,
                Utils.hasExistingInput(meta, Constants.INPUT.IDX_DNA_TUMOR) ? [Utils.getInput(meta, Constants.INPUT.IDX_DNA_TUMOR)] : idxs,
            ]
        }
        .branch { meta, alns, idxs ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.REDUX_DIR_TUMOR)
            runnable: alns && ! has_existing
            skip: true
                return meta
        }

    ch_inputs_normal_sorted = ch_dna_normal
        .map { meta, alns, idxs ->
            return [
                meta,
                Utils.hasExistingInput(meta, Constants.INPUT.ALN_DNA_NORMAL) ? [Utils.getInput(meta, Constants.INPUT.ALN_DNA_NORMAL)] : alns,
                Utils.hasExistingInput(meta, Constants.INPUT.IDX_DNA_NORMAL) ? [Utils.getInput(meta, Constants.INPUT.IDX_DNA_NORMAL)] : idxs,
            ]
        }
        .branch { meta, alns, idxs ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.REDUX_DIR_NORMAL)
            runnable: alns && ! has_existing
            skip: true
                return meta
        }

    ch_inputs_donor_sorted = ch_dna_donor
        .map { meta, alns, idxs ->
            return [
                meta,
                Utils.hasExistingInput(meta, Constants.INPUT.ALN_DNA_DONOR) ? [Utils.getInput(meta, Constants.INPUT.ALN_DNA_DONOR)] : alns,
                Utils.hasExistingInput(meta, Constants.INPUT.IDX_DNA_DONOR) ? [Utils.getInput(meta, Constants.INPUT.IDX_DNA_DONOR)] : idxs,
            ]
        }
        .branch { meta, alns, idxs ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.REDUX_DIR_DONOR)
            runnable: alns && ! has_existing
            skip: true
            return meta
        }

    // Create process input channel
    // channel: [ meta_redux, [aln, ...], [idx, ...] ]
    ch_redux_inputs = channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable.map { meta, alns, idxs -> [meta, Utils.getTumorDnaSample(meta), 'tumor', alns, idxs] },
            ch_inputs_normal_sorted.runnable.map { meta, alns, idxs -> [meta, Utils.getNormalDnaSample(meta), 'normal', alns, idxs] },
            ch_inputs_donor_sorted.runnable.map { meta, alns, idxs -> [meta, Utils.getDonorDnaSample(meta), 'donor', alns, idxs] },
        )
        .multiMap { meta, meta_sample, sample_type, alns, idxs ->

            def sample_id = meta_sample.getOrDefault('longitudinal_sample_id', meta_sample['sample_id'])

            def meta_redux = [
                key: meta.group_id,
                id: "${meta.group_id}_${sample_id}",
                sample_id: sample_id,
                sample_type: sample_type,
            ]

            sample_data: [meta_redux, alns, idxs]
            generate_tsvs_only: meta_sample.getOrDefault(Constants.InfoField.GENERATE_REDUX_TSVS_ONLY, false)
        }

    // Run process
    REDUX(
        ch_redux_inputs.sample_data,
        genome_fasta,
        genome_version,
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
    )

    // Split into sample type channels
    // channel: [ meta_redux, redux_dir ]
    ch_redux_out_sorted = channel.topic('redux_dir')
        .branch { meta_redux, redux_dir ->
            assert ['tumor', 'normal', 'donor'].contains(meta_redux.sample_type)
            tumor: meta_redux.sample_type == 'tumor'
            normal: meta_redux.sample_type == 'normal'
            donor: meta_redux.sample_type == 'donor'
            placeholder: true
        }

    // Set outputs, restoring original meta
    // channel: [ meta, redux_dir ]
    ch_outputs_tumor = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_redux_out_sorted.tumor, ch_inputs),
            ch_inputs_tumor_sorted.skip.map { meta -> [meta, []] },
        )

    // channel: [ meta, redux_dir ]
    ch_outputs_normal = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_redux_out_sorted.normal, ch_inputs),
            ch_inputs_normal_sorted.skip.map { meta -> [meta, []] },
        )

    // channel: [ meta, redux_dir ]
    ch_outputs_donor = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_redux_out_sorted.donor, ch_inputs),
            ch_inputs_donor_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    tumor_dir  = ch_outputs_tumor  // channel: [ meta, redux_dir ]
    normal_dir = ch_outputs_normal // channel: [ meta, redux_dir ]
    donor_dir  = ch_outputs_donor  // channel: [ meta, redux_dir ]
}
