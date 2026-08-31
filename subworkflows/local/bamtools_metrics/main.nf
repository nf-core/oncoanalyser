//
// Bam Tools calculates summary statistics for BAMs
//

nextflow.enable.types = true

include { BAMTOOLS } from '../../../modules/local/bamtools/main'

include { getNormalReduxDirAlignment } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { getTumorReduxDirAlignment  } from '../utils_nfcore_oncoanalyser_pipeline/accessors_alignments'
include { getInput                   } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalDnaSample         } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSample          } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasInput                   } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { groupByMeta                } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { joinMeta                   } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { restoreMeta                } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { FileType                   } from '../utils_nfcore_oncoanalyser_pipeline/types_enums'
include { selectCurrentOrExisting    } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow BAMTOOLS_METRICS {
    take:
    // Sample data
    ch_inputs             : Channel<Map>              // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor    : Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal   : Channel<Tuple<Map, Path>> // channel: [mandatory] [ meta, redux_dir ]

    // Reference data
    genome_fasta          : Channel<Path>             // channel: [mandatory] /path/to/genome_fasta
    genome_version        : Channel<String>           // channel: [mandatory] genome version
    genome_fai            : Channel<Path>             // channel: [mandatory] /path/to/genome_fai
    driver_gene_panel     : Channel<Path>             // channel: [mandatory] /path/to/driver_gene_panel
    ensembl_data_resources: Channel<Path>             // channel: [mandatory] /path/to/ensembl_data_resources/
    target_regions_bed    : Channel<Path>?            // channel: [optional]  /path/to/target_regions_bed

    main:
    // Select input sources then sort
    // channel: runnable: [ meta, aln, idx ]
    // channel: skip: [ meta ]
    ch_inputs_tumor_sorted = ch_redux_dir_tumor
        .map { meta, redux_dir ->

            def redux_dir_selected = selectCurrentOrExisting(redux_dir, getInput(getTumorDnaSample(meta), FileType.REDUX_DIR))
            def (aln, idx) = getTumorReduxDirAlignment(meta, redux_dir_selected)

            return [meta, aln, idx]

        }
        .branch { meta, aln, idx ->
            def has_existing = hasInput(getTumorDnaSample(meta), FileType.BAMTOOLS_DIR)
            runnable: aln && ! has_existing
            skip: true
                return meta
        }

    // channel: runnable: [ meta, aln, idx ]
    // channel: skip: [ meta ]
    ch_inputs_normal_sorted = ch_redux_dir_normal
        .map { meta, redux_dir ->

            def redux_dir_selected = selectCurrentOrExisting(redux_dir, getInput(getNormalDnaSample(meta), FileType.REDUX_DIR))
            def (aln, idx) = getNormalReduxDirAlignment(meta, redux_dir_selected)

            return [meta, aln, idx]
        }
        .branch { meta, aln, idx ->
            def has_existing = hasInput(getNormalDnaSample(meta), FileType.BAMTOOLS_DIR)
            runnable: aln && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_bamtools, aln, idx ]
    ch_bamtools_inputs = channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable.map { meta, aln, idx -> [meta, getTumorDnaSample(meta), 'tumor', aln, idx] },
            ch_inputs_normal_sorted.runnable.map { meta, aln, idx -> [meta, getNormalDnaSample(meta), 'normal', aln, idx] },
        )
        .map { meta, meta_sample, sample_type, aln, idx ->

            def meta_bamtools = record(
                key: meta.case_id,
                id: "${meta.case_id}:${meta_sample.sample_id}",
                sample_id: meta_sample.sample_id,
                sample_type: sample_type,
            )

            return [meta_bamtools, aln, idx]
        }

    // Run process
    BAMTOOLS(
        ch_bamtools_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        driver_gene_panel,
        ensembl_data_resources,
        target_regions_bed,
    )

    // Sort into a tumor and normal channel
    ch_bamtools_out = channel.topic('bamtools_metrics_dir')
        .branch { meta_bamtools, bamtools_dir ->
            assert ['tumor', 'normal'].contains(meta_bamtools.sample_type)
            tumor: meta_bamtools.sample_type == 'tumor'
            normal: meta_bamtools.sample_type == 'normal'
            placeholder: true
        }

    // Set outputs, restoring original meta
    // channel: [ meta, bamtools_dir ]
    ch_tumor_out = channel.empty()
        .mix(
            restoreMeta(ch_bamtools_out.tumor, ch_inputs),
            ch_inputs_tumor_sorted.skip.map { meta -> [meta, null] },
        )

    // channel: [ meta, bamtools_dir ]
    ch_normal_out = channel.empty()
        .mix(
            restoreMeta(ch_bamtools_out.normal, ch_inputs),
            ch_inputs_normal_sorted.skip.map { meta -> [meta, null] },
        )

    emit:
    tumor_dir  = ch_tumor_out  // channel: [ meta, bamtools_dir ]
    normal_dir = ch_normal_out // channel: [ meta, bamtools_dir ]
}
