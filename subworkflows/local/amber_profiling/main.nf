//
// AMBER determines b-allele frequencies at predetermined positions
//

nextflow.enable.types = true

include { AMBER  } from '../../../modules/local/amber/main'

include { FileType                    } from '../utils_nfcore_oncoanalyser_pipeline/types_enums'
include { groupByMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { joinMeta                    } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { restoreMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { getDonorDnaSample           } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getDonorDnaSampleNames      } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getDonorReduxDirAlignments  } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { getInput                    } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getLongitudinalSampleName } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalDnaSample          } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalDnaSampleName      } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalReduxDirAlignment  } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { getTumorDnaSample           } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSampleName       } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorReduxDirAlignment   } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { hasInput                    } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { selectCurrentOrExisting     } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow AMBER_PROFILING {
    take:
    // Sample data
    ch_inputs: Channel<Map>            // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor: Channel<Tuple<Map, Path>>   // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal: Channel<Tuple<Map, Path>>  // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_donor: Channel<Tuple<Map, Path>>   // channel: [mandatory] [ meta, redux_dir ]

    // Reference data
    genome_fasta: Channel<Path>         // channel: [mandatory] /path/to/genome_fasta
    genome_version: Channel<String>       // channel: [mandatory] genome version
    genome_fai: Channel<Path>           // channel: [mandatory] /path/to/genome_fai
    heterozygous_sites: Channel<Path>   // channel: [optional]  /path/to/heterozygous_sites
    target_regions_bed: Channel<Path>   // channel: [optional]  /path/to/target_regions_bed
    tumor_min_depth: Integer?      // integer: [optional]  -tumor_min_depth argument value

    // Params
    sequencing_platform: String  // string:  [mandatory] sequencing platform
    purity_estimate_mode: Boolean // boolean: [mandatory] Set purity estimate mode

    main:
    // Select input sources then sort
    // channel: runnable: [ meta, tumor_aln, tumor_idx, normal_aln, normal_idx ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = groupByMeta([
        ch_redux_dir_tumor,
        ch_redux_dir_normal,
        ch_redux_dir_donor,
    ])
        .map { meta, redux_dir_tumor, redux_dir_normal, redux_dir_donor ->

            def redux_dir_tumor_selected = selectCurrentOrExisting(redux_dir_tumor, getInput(getTumorDnaSample(meta), FileType.REDUX_DIR))
            def redux_dir_normal_selected = selectCurrentOrExisting(redux_dir_normal, getInput(getNormalDnaSample(meta), FileType.REDUX_DIR))
            def redux_dir_donor_selected = selectCurrentOrExisting(redux_dir_donor, getInput(getDonorDnaSample(meta), FileType.REDUX_DIR))

            def (tumor_aln, tumor_idx) = getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)
            def (normal_aln, normal_idx) = getNormalReduxDirAlignment(meta, redux_dir_normal_selected)
            def donor_alignments = getDonorReduxDirAlignments(meta, redux_dir_donor_selected)
            def donor_alns = donor_alignments.collect { aln, idx -> aln }
            def donor_idxs = donor_alignments.collect { aln, idx -> idx }

            return [meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_alns, donor_idxs]

        }
        .branch { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_alns, donor_idxs ->

            def has_existing = hasInput(getTumorDnaSample(meta), FileType.AMBER_DIR)
            def runnable_standard = ! purity_estimate_mode && tumor_aln && ! has_existing

            // TODO(SW): must improve handling through separation of sample information in meta; currently unable to provide ccfDNA AMBER directory in samplesheet
            def runnable_purity_estimate = purity_estimate_mode && normal_aln

            runnable: runnable_standard || runnable_purity_estimate
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_amber, tumor_aln, normal_aln, donor_aln, tumor_idx, normal_idx, donor_idx ]
    ch_amber_inputs = ch_inputs_sorted.runnable
        .map { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_alns, donor_idxs ->

            def tumor_id
            if (purity_estimate_mode) {
                tumor_id = getLongitudinalSampleName(meta)
            } else {
                tumor_id = getTumorDnaSampleName(meta)
            }

            def meta_amber = record(
                key: meta.case_id,
                id: meta.case_id,
                tumor_id: tumor_id,
                normal_id: normal_aln ? getNormalDnaSampleName(meta) : null,
                donor_ids: donor_alns ? getDonorDnaSampleNames(meta) : null,
            )

            return [meta_amber, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_alns, donor_idxs]
        }

    // Run process
    AMBER(
        ch_amber_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        heterozygous_sites,
        target_regions_bed,
        tumor_min_depth,
        sequencing_platform,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, amber_dir ]
    ch_outputs = channel.empty()
        .mix(
            restoreMeta(channel.topic('amber_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    emit:
    amber_dir = ch_outputs // channel: [ meta, amber_dir ]
}
