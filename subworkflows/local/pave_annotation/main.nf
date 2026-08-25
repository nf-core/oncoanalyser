//
// PAVE annotates somatic and germline variant VCFs with gene and transcript coding and protein effects
//

nextflow.enable.types = true

include { PAVE_GERMLINE } from '../../../modules/local/pave/germline/main'
include { PAVE_SOMATIC  } from '../../../modules/local/pave/somatic/main'

include { getSageGermlineVcf      } from '../utils_nfcore_oncoanalyser_pipeline/accessors_outputs'
include { getSageSomaticVcf       } from '../utils_nfcore_oncoanalyser_pipeline/accessors_outputs'
include { getVcfTbi               } from '../utils_nfcore_oncoanalyser_pipeline/accessors_outputs'
include { getInput                } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalDnaSample      } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSample       } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSampleName   } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasInput                } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasNormalDna            } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { hasTumorDna             } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { FileType                } from '../utils_nfcore_oncoanalyser_pipeline/types_enums'
include { groupByMeta             } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { joinMeta                } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { restoreMeta             } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { selectCurrentOrExisting } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow PAVE_ANNOTATION {
    take:
    // Sample data
    ch_inputs: Channel<Map>              // channel: [mandatory] [ meta ]
    ch_sage_dir_somatic: Channel<Tuple<Map, Path>>    // channel: [mandatory] [ meta, sage_dir ]
    ch_sage_dir_germline: Channel<Tuple<Map, Path>>   // channel: [mandatory] [ meta, sage_dir ]

    // Reference data
    genome_fasta: Channel<Path>           // channel: [mandatory] /path/to/genome_fasta
    genome_version: Channel<String>         // channel: [mandatory] genome version
    genome_fai: Channel<Path>             // channel: [mandatory] /path/to/genome_fai
    pon_artefacts: Channel<Path>          // channel: [optional]  /path/to/pon_artefacts
    sage_pon: Channel<Path>               // channel: [mandatory] /path/to/sage_pon
    sage_blocklist_regions: Channel<Path> // channel: [mandatory] /path/to/sage_blocklist_regions
    sage_blocklist_sites: Channel<Path>   // channel: [mandatory] /path/to/sage_blocklist_sites
    clinvar_annotations: Channel<Path>    // channel: [mandatory] /path/to/clinvar_annotations
    segment_mappability: Channel<Path>    // channel: [mandatory] /path/to/segment_mappability
    driver_gene_panel: Channel<Path>      // channel: [mandatory] /path/to/driver_gene_panel
    ensembl_data_resources: Channel<Path> // channel: [mandatory] /path/to/ensembl_data_resources/
    gnomad_resource: Channel<Path>        // channel: [mandatory] /path/to/gnomad_resource

    // Params
    sequencing_platform: String    // string:  [mandatory] sequencing platform

    main:
    //
    // MODULE: PAVE germline
    //
    // Select input sources then sort
    // channel: runnable: [ meta, sage_vcf, sage_tbi ]
    // channel: skip: [ meta ]
    ch_sage_germline_inputs_sorted = ch_sage_dir_germline
        .map { meta, sage_dir ->

            def sage_dir_selected = selectCurrentOrExisting(sage_dir, getInput(getNormalDnaSample(meta), FileType.SAGE_DIR))
            def sage_vcf = getSageGermlineVcf(getTumorDnaSampleName(meta), sage_dir_selected)
            def sage_tbi = getVcfTbi(sage_vcf)

            return [meta, sage_vcf, sage_tbi]
        }
        .branch { meta, sage_vcf, sage_tbi ->

            def has_existing = hasInput(getNormalDnaSample(meta), FileType.PAVE_DIR)

            runnable: hasTumorDna(meta) && hasNormalDna(meta) && sage_vcf && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_pave, sage_vcf, sage_tbi ]
    ch_pave_germline_inputs = ch_sage_germline_inputs_sorted.runnable
        .map { meta, sage_vcf, sage_tbi ->

            def meta_pave = record(
                key: meta.case_id,
                id: meta.case_id,
                sample_id: getTumorDnaSampleName(meta),
            )

            return [meta_pave, sage_vcf, sage_tbi]
        }

    // Run process
    PAVE_GERMLINE(
        ch_pave_germline_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        sage_blocklist_regions,
        sage_blocklist_sites,
        clinvar_annotations,
        segment_mappability,
        driver_gene_panel,
        ensembl_data_resources,
        sequencing_platform,
    )

    //
    // MODULE: PAVE somatic
    //
    // Select input sources then sort
    // channel: runnable: [ meta, sage_vcf, sage_tbi ]
    // channel: skip: [ meta ]
    ch_sage_somatic_inputs_sorted = ch_sage_dir_somatic
        .map { meta, sage_dir ->

            def sage_dir_selected = selectCurrentOrExisting(sage_dir, getInput(getTumorDnaSample(meta), FileType.SAGE_DIR))
            def sage_vcf = getSageSomaticVcf(getTumorDnaSampleName(meta), sage_dir_selected)
            def sage_tbi = getVcfTbi(sage_vcf)

            return [meta, sage_vcf, sage_tbi]
        }
        .branch { meta, sage_vcf, sage_tbi ->

            def has_existing = hasInput(getTumorDnaSample(meta), FileType.PAVE_DIR)

            runnable: hasTumorDna(meta) && sage_vcf && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_pave, sage_vcf, sage_tbi ]
    ch_pave_somatic_inputs = ch_sage_somatic_inputs_sorted.runnable
        .map { meta, sage_vcf, sage_tbi ->

            def meta_pave = record(
                key: meta.case_id,
                id: meta.case_id,
                sample_id: getTumorDnaSampleName(meta),
            )

            return [meta_pave, sage_vcf, sage_tbi]
        }

    // Run process
    PAVE_SOMATIC(
        ch_pave_somatic_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        pon_artefacts,
        sage_pon,
        clinvar_annotations,
        segment_mappability,
        driver_gene_panel,
        ensembl_data_resources,
        gnomad_resource,
        sequencing_platform,
    )

    //
    // STEP: Handle outputs
    //
    // Set outputs, restoring original meta
    // channel: [ meta, pave_dir ]
    ch_outputs_somatic = channel.empty()
        .mix(
            restoreMeta(channel.topic('pave_somatic_dir'), ch_inputs),
            ch_sage_somatic_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    // channel: [ meta, pave_dir ]
    ch_outputs_germline = channel.empty()
        .mix(
            restoreMeta(channel.topic('pave_germline_dir'), ch_inputs),
            ch_sage_germline_inputs_sorted.skip.map { meta -> [meta, null] },
        )

    emit:
    germline_dir = ch_outputs_germline // channel: [ meta, pave_dir ]
    somatic_dir  = ch_outputs_somatic  // channel: [ meta, pave_dir ]
}
