//
// PURPLE is a CNV caller that infers purity/ploidy and recovers low-confidence SVs
//

nextflow.enable.types = true

include { PURPLE  } from '../../../modules/local/purple/main'

include { FileType                 } from '../utils_nfcore_oncoanalyser_pipeline/types'
include { groupByMeta              } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { joinMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { restoreMeta              } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { getInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getNormalDnaSample       } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getNormalDnaSampleName   } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSample        } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSampleName    } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorReduxTsvs        } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { hasInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { hasNormalDna             } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { selectCurrentOrExisting  } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow PURPLE_CALLING {
    take:
    // Sample data
    ch_inputs: Channel<Map>                    // channel: [mandatory] [ meta ]
    ch_amber_dir: Channel<Tuple<Map, Path>>                 // channel: [mandatory] [ meta, amber_dir ]
    ch_cobalt_dir: Channel<Tuple<Map, Path>>                // channel: [mandatory] [ meta, cobalt_dir ]
    ch_esvee_dir: Channel<Tuple<Map, Path>>                 // channel: [mandatory] [ meta, esvee_dir ]
    ch_pave_somatic_dir: Channel<Tuple<Map, Path>>          // channel: [mandatory] [ meta, pave_dir ]
    ch_pave_germline_dir: Channel<Tuple<Map, Path>>         // channel: [mandatory] [ meta, pave_dir ]
    ch_redux_dir_tumor: Channel<Tuple<Map, Path>>           // channel: [optional]  [ meta, redux_dir ]

    // Reference data
    genome_fasta: Channel<Path>                 // channel: [mandatory] /path/to/genome_fasta
    genome_version: Channel<String>               // channel: [mandatory] genome version
    genome_fai: Channel<Path>                   // channel: [mandatory] /path/to/genome_fai
    genome_dict: Channel<Path>                  // channel: [mandatory] /path/to/genome_dict
    gc_profile: Channel<Path>                   // channel: [mandatory] /path/to/gc_profile
    sage_known_hotspots_somatic: Channel<Path>  // channel: [mandatory] /path/to/sage_known_hotspots_somatic
    sage_known_hotspots_germline: Channel<Path> // channel: [optional]  /path/to/sage_known_hotspots_germline
    driver_gene_panel: Channel<Path>            // channel: [mandatory] /path/to/driver_gene_panel
    ensembl_data_resources: Channel<Path>       // channel: [mandatory] /path/to/ensembl_data_resources/
    germline_amp_del_freq: Channel<Path>        // channel: [optional]  /path/to/germline_amp_del_freq
    target_regions_bed: Channel<Path>           // channel: [optional]  /path/to/target_regions_bed

    main:
    // Select input sources then sort
    // channel: runnable: [ meta, amber_dir, cobalt_dir, esvee_dir, pave_somatic_dir, pave_germline_dir, [redux_tsv_tumor, ...] ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = groupByMeta([
        ch_amber_dir,
        ch_cobalt_dir,
        ch_esvee_dir,
        ch_pave_somatic_dir,
        ch_pave_germline_dir,
        ch_redux_dir_tumor,
    ])
        .map { meta, amber_dir, cobalt_dir, esvee_dir, pave_somatic_dir, pave_germline_dir, redux_dir_tumor ->

            def tumor_dir_selected = selectCurrentOrExisting(redux_dir_tumor, getInput(getTumorDnaSample(meta), FileType.REDUX_DIR))
            def tumor_tsvs = getTumorReduxTsvs(meta, tumor_dir_selected)

            return [
                meta,
                selectCurrentOrExisting(amber_dir, getInput(getTumorDnaSample(meta), FileType.AMBER_DIR)),
                selectCurrentOrExisting(cobalt_dir, getInput(getTumorDnaSample(meta), FileType.COBALT_DIR)),
                selectCurrentOrExisting(esvee_dir, getInput(getTumorDnaSample(meta), FileType.ESVEE_DIR)),
                selectCurrentOrExisting(pave_somatic_dir, getInput(getTumorDnaSample(meta), FileType.PAVE_DIR)),
                selectCurrentOrExisting(pave_germline_dir, getInput(getNormalDnaSample(meta), FileType.PAVE_DIR)),
                tumor_tsvs.findAll { f -> f.exists() },
            ]
        }
        .branch { meta, amber_dir, cobalt_dir, esvee_dir, pave_somatic_dir, pave_germline_dir, redux_tsvs_tumor ->

            def has_existing = hasInput(getTumorDnaSample(meta), FileType.PURPLE_DIR)

            runnable: amber_dir && cobalt_dir && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_purple, amber_dir, cobalt_dir, esvee_dir, pave_somatic_dir, pave_germline_dir, [redux_tsv_tumor, ...] ]
    ch_purple_inputs = ch_inputs_sorted.runnable
        .map { meta, amber_dir, cobalt_dir, esvee_dir, pave_somatic_dir, pave_germline_dir, redux_tsvs_tumor ->

            def meta_purple = record(
                key: meta.case_id,
                id: meta.case_id,
                tumor_id: getTumorDnaSampleName(meta),
                normal_id: hasNormalDna(meta) ? getNormalDnaSampleName(meta) : null,
            )

            return [meta_purple, amber_dir, cobalt_dir, esvee_dir, pave_somatic_dir, pave_germline_dir, redux_tsvs_tumor]
        }

    // Run process
    PURPLE(
        ch_purple_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        gc_profile,
        sage_known_hotspots_somatic,
        sage_known_hotspots_germline,
        driver_gene_panel,
        ensembl_data_resources,
        germline_amp_del_freq,
        target_regions_bed,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, purple_dir ]
    ch_outputs = channel.empty()
        .mix(
            restoreMeta(channel.topic('purple_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    emit:
    purple_dir = ch_outputs // channel: [ meta, purple_dir ]
}
