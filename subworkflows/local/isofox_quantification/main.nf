//
// Isofox estimates transcript abundance, detects novel SJs, and identifies fusion events
//

nextflow.enable.types = true

include { ISOFOX  } from '../../../modules/local/isofox/run/main'

include { FileType                 } from '../utils_nfcore_oncoanalyser_pipeline/types_enums'
include { groupByMeta              } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { joinMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { restoreMeta              } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { getInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSampleName    } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorRnaSample        } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorRnaSampleName    } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { selectCurrentOrExisting  } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow ISOFOX_QUANTIFICATION {
    take:
    // Sample data
    ch_inputs: Channel<Map>                  // channel: [mandatory] [ meta ]
    ch_tumor_rna_aln: Channel<Tuple<Map, Path, Path>>           // channel: [mandatory] [ meta, aln, idx ]

    // Reference data
    genome_fasta: Channel<Path>               // channel: [mandatory] /path/to/genome_fasta
    genome_version: Channel<String>             // channel: [mandatory] genome version
    genome_fai: Channel<Path>                 // channel: [mandatory] /path/to/genome_fai
    ensembl_data_resources: Channel<Path>     // channel: [mandatory] /path/to/ensembl_data_resources/
    driver_gene_panel: Channel<Path>          // channel: [mandatory] /path/to/driver_gene_panel
    known_fusion_data: Channel<Path>          // channel: [mandatory] /path/to/known_fusion_data
    isofox_excluded_regions: Channel<Path>    // channel: [mandatory] /path/to/isofox_excluded_regions
    isofox_gene_distribution: Channel<Path>   // channel: [mandatory] /path/to/isofox_gene_distribution
    isofox_alt_sj_distribution: Channel<Path> // channel: [mandatory] /path/to/isofox_alt_sj_distribution
    isofox_counts: Channel<Path>              // channel: [mandatory] /path/to/isofox_counts
    isofox_gc_ratios: Channel<Path>           // channel: [mandatory] /path/to/isofox_gc_ratios
    isofox_tpm_norm: Channel<Path>            // channel: [optional]  /path/to/isofox_tpm_norm

    // Params
    isofox_functions: String?       //  string: [optional]  Isofox functions
    isofox_read_length: Integer     //  string: [mandatory] Isofox read length

    main:
    // Select input sources then sort
    // channel: runnable: [ meta, tumor_aln, tumor_idx ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_tumor_rna_aln
        .map { meta, tumor_aln, tumor_idx ->
            return [
                meta,
                selectCurrentOrExisting(tumor_aln, getInput(getTumorRnaSample(meta), FileType.ALN)),
                selectCurrentOrExisting(tumor_idx, getInput(getTumorRnaSample(meta), FileType.IDX)),
            ]
        }
        .branch { meta, tumor_aln, tumor_idx ->
            def has_existing = hasInput(getTumorRnaSample(meta), FileType.ISOFOX_DIR)
            runnable: tumor_aln && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_isofox, tumor_aln, tumor_idx ]
    ch_isofox_inputs = ch_inputs_sorted.runnable
        .map { meta, tumor_aln, tumor_idx ->

            def meta_isofox = record(
                key: meta.case_id,
                id: meta.case_id,
                sample_id: getTumorDnaSampleName(meta) ?: getTumorRnaSampleName(meta),
            )

            return [meta_isofox, tumor_aln, tumor_idx]
        }

    // Run process
    ISOFOX(
        ch_isofox_inputs,
        isofox_functions,
        isofox_read_length,
        genome_fasta,
        genome_version,
        genome_fai,
        ensembl_data_resources,
        driver_gene_panel,
        known_fusion_data,
        isofox_excluded_regions,
        isofox_gene_distribution,
        isofox_alt_sj_distribution,
        isofox_counts,
        isofox_gc_ratios,
        isofox_tpm_norm,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, isofox_dir ]
    ch_outputs = channel.empty()
        .mix(
            restoreMeta(channel.topic('isofox_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    emit:
    isofox_dir = ch_outputs // channel: [ meta, isofox_dir ]
}
