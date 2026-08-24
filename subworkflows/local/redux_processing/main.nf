//
// Apply post-alignment processing
//

nextflow.enable.types = true

include { REDUX  } from '../../../modules/local/redux/main'

include { FileType            } from '../utils_nfcore_oncoanalyser_pipeline/types'
include { groupByMeta         } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { joinMeta            } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { restoreMeta         } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { getDonorDnaSample   } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getDonorDnaSamples  } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getInput            } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getNormalDnaAln     } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getNormalDnaSample  } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaAln      } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSample   } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { hasInput            } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { hasNormalDnaAln     } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { hasTumorDnaAln      } from '../utils_nfcore_oncoanalyser_pipeline/accessors'

workflow REDUX_PROCESSING {
    take:
    // Sample data
    ch_inputs: Channel<Map>              // channel: [mandatory] [ meta ]
    ch_dna_tumor: Channel<Tuple<Map, List<Path>, List<Path>>>           // channel: [mandatory] [ meta, [aln, ...], [idx, ...] ]
    ch_dna_normal: Channel<Tuple<Map, List<Path>, List<Path>>>          // channel: [mandatory] [ meta, [aln, ...], [idx, ...] ]
    ch_dna_donor: Channel<Tuple<Map, List<Path>, List<Path>>>           // channel: [mandatory] [ meta, [aln, ...], [idx, ...] ]

    // Reference data
    genome_fasta: Channel<Path>           // channel: [mandatory] /path/to/genome_fasta
    genome_version: Channel<String>         // channel: [mandatory] genome version
    genome_fai: Channel<Path>             // channel: [mandatory] /path/to/genome_fai
    genome_dict: Channel<Path>            // channel: [mandatory] /path/to/genome_dict
    unmap_regions: Channel<Path>          // channel: [mandatory] /path/to/unmap_regions
    msi_jitter_sites: Channel<Path>       // channel: [mandatory] /path/to/msi_jitter_sites
    msi_model_coefficients: Channel<Path> // channel: [mandatory] /path/to/msi_model_coefficients
    msi_model_error_rates: Channel<Path>  // channel: [mandatory] /path/to/msi_model_error_rates

    // Params
    sequencing_platform: String    // string:  [mandatory] sequencing platform
    targeted_mode: Boolean          // boolean: [mandatory] Set targeted mode
    umi_enable: Boolean             // boolean: [mandatory] enable UMI processing
    umi_duplex_delim: String?       // string:  [optional] UMI duplex delimiter

    main:
    // Select input sources then sort, separating by sample type
    // channel: runnable: [ meta, [aln, ...], [idx, ...] ]
    // channel: skip: [ meta ]
    ch_inputs_tumor_sorted = ch_dna_tumor
        .map { meta, alns, idxs ->
            return [
                meta,
                hasTumorDnaAln(meta) ? [getTumorDnaAln(meta)] : alns,
                hasInput(getTumorDnaSample(meta), FileType.IDX) ? [getInput(getTumorDnaSample(meta), FileType.IDX)] : idxs,
            ]
        }
        .branch { meta, alns, idxs ->
            def has_existing = hasInput(getTumorDnaSample(meta), FileType.REDUX_DIR)
            runnable: alns && ! has_existing
            skip: true
                return meta
        }

    ch_inputs_normal_sorted = ch_dna_normal
        .map { meta, alns, idxs ->
            return [
                meta,
                hasNormalDnaAln(meta) ? [getNormalDnaAln(meta)] : alns,
                hasInput(getNormalDnaSample(meta), FileType.IDX) ? [getInput(getNormalDnaSample(meta), FileType.IDX)] : idxs,
            ]
        }
        .branch { meta, alns, idxs ->
            def has_existing = hasInput(getNormalDnaSample(meta), FileType.REDUX_DIR)
            runnable: alns && ! has_existing
            skip: true
                return meta
        }

    ch_inputs_donor_sorted = ch_dna_donor
        .map { meta, sample_id, alns, idxs ->
            def donor = sample_id ? getDonorDnaSamples(meta).find { it.sample_id == sample_id } : getDonorDnaSample(meta)
            def donor_aln = donor ? getInput(donor, FileType.ALN) ?: getInput(donor, FileType.ALN_REDUX) : null
            return [
                meta,
                donor,
                donor_aln ? [donor_aln] : alns,
                donor && hasInput(donor, FileType.IDX) ? [getInput(donor, FileType.IDX)] : idxs,
            ]
        }
        .branch { meta, donor, alns, idxs ->
            def has_existing = donor && hasInput(donor, FileType.REDUX_DIR)
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
            ch_inputs_donor_sorted.runnable.map { meta, donor, alns, idxs -> [meta, donor, 'donor', alns, idxs] },
        )
        .multiMap { meta, meta_sample, sample_type, alns, idxs ->

            def sample_id = meta_sample.sample_id

            def meta_redux = record(
                key: meta.case_id,
                id: "${meta.case_id}_${sample_id}",
                sample_id: sample_id,
                sample_type: sample_type,
            )

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
            ch_inputs_tumor_sorted.skip.map { meta -> [meta, null] },
        )

    // channel: [ meta, redux_dir ]
    ch_outputs_normal = channel.empty()
        .mix(
            restoreMeta(ch_redux_out_sorted.normal, ch_inputs),
            ch_inputs_normal_sorted.skip.map { meta -> [meta, null] },
        )

    // channel: [ meta, [redux_dir, ...] ]
    ch_outputs_donor = channel.empty()
        .mix(
            restoreMeta(ch_redux_out_sorted.donor, ch_inputs).groupTuple(),
            ch_inputs_donor_sorted.skip.map { meta -> [meta, null] },
        )

    emit:
    tumor_dir  = ch_outputs_tumor  // channel: [ meta, redux_dir ]
    normal_dir = ch_outputs_normal // channel: [ meta, redux_dir ]
    donor_dir  = ch_outputs_donor  // channel: [ meta, [redux_dir, ...] ]
}
