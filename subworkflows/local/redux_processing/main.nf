//
// Apply post-alignment processing
//

include { REDUX } from '../../../modules/local/redux/main'

include { groupByMeta          } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { joinMeta             } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { restoreMeta          } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { getDonorDnaAln       } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getDonorDnaBai       } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getDonorDnaSample    } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getNormalDnaAln      } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getNormalDnaBai      } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getNormalDnaSample   } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getTumorDnaAln       } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getTumorDnaBai       } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { getTumorDnaSample    } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { hasDonorDnaAln       } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { hasDonorDnaBai       } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { hasDonorDnaReduxDir  } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { hasNormalDnaAln      } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { hasNormalDnaBai      } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { hasNormalDnaReduxDir } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { hasTumorDnaAln       } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { hasTumorDnaBai       } from '../utils_nfcore_oncoanalyser_pipeline/utils'
include { hasTumorDnaReduxDir  } from '../utils_nfcore_oncoanalyser_pipeline/utils'

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
                hasTumorDnaAln(meta) ? [getTumorDnaAln(meta)] : alns,
                hasTumorDnaBai(meta) ? [getTumorDnaBai(meta)] : idxs,
            ]
        }
        .branch { meta, alns, idxs ->
            def has_existing = hasTumorDnaReduxDir(meta)
            runnable: alns && ! has_existing
            skip: true
                return meta
        }

    ch_inputs_normal_sorted = ch_dna_normal
        .map { meta, alns, idxs ->
            return [
                meta,
                hasNormalDnaAln(meta) ? [getNormalDnaAln(meta)] : alns,
                hasNormalDnaBai(meta) ? [getNormalDnaBai(meta)] : idxs,
            ]
        }
        .branch { meta, alns, idxs ->
            def has_existing = hasNormalDnaReduxDir(meta)
            runnable: alns && ! has_existing
            skip: true
                return meta
        }

    ch_inputs_donor_sorted = ch_dna_donor
        .map { meta, alns, idxs ->
            return [
                meta,
                hasDonorDnaAln(meta) ? [getDonorDnaAln(meta)] : alns,
                hasDonorDnaBai(meta) ? [getDonorDnaBai(meta)] : idxs,
            ]
        }
        .branch { meta, alns, idxs ->
            def has_existing = hasDonorDnaReduxDir(meta)
            runnable: alns && ! has_existing
            skip: true
            return meta
        }

    // Create process input channel
    // channel: [ meta_redux, [aln, ...], [idx, ...] ]
    ch_redux_inputs = channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable.map { meta, alns, idxs -> [meta, getTumorDnaSample(meta), 'tumor', alns, idxs] },
            ch_inputs_normal_sorted.runnable.map { meta, alns, idxs -> [meta, getNormalDnaSample(meta), 'normal', alns, idxs] },
            ch_inputs_donor_sorted.runnable.map { meta, alns, idxs -> [meta, getDonorDnaSample(meta), 'donor', alns, idxs] },
        )
        .multiMap { meta, meta_sample, sample_type, alns, idxs ->

            def sample_id = meta_sample.sample_id

            def meta_redux = [
                key: meta.case_id,
                id: "${meta.case_id}_${sample_id}",
                sample_id: sample_id,
                sample_type: sample_type,
            ]

            sample_data: [meta_redux, alns, idxs]
            generate_tsvs_only: meta_sample.generate_redux_tsvs_only
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
            restoreMeta(ch_redux_out_sorted.tumor, ch_inputs),
            ch_inputs_tumor_sorted.skip.map { meta -> [meta, []] },
        )

    // channel: [ meta, redux_dir ]
    ch_outputs_normal = channel.empty()
        .mix(
            restoreMeta(ch_redux_out_sorted.normal, ch_inputs),
            ch_inputs_normal_sorted.skip.map { meta -> [meta, []] },
        )

    // channel: [ meta, redux_dir ]
    ch_outputs_donor = channel.empty()
        .mix(
            restoreMeta(ch_redux_out_sorted.donor, ch_inputs),
            ch_inputs_donor_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    tumor_dir  = ch_outputs_tumor  // channel: [ meta, redux_dir ]
    normal_dir = ch_outputs_normal // channel: [ meta, redux_dir ]
    donor_dir  = ch_outputs_donor  // channel: [ meta, redux_dir ]
}
