//
// PAVE annotates somatic and germline variant VCFs with gene and transcript coding and protein effects
//

nextflow.enable.types = true

include { PAVE_GERMLINE  } from '../../../modules/local/pave/germline/main'
include { PAVE_SOMATIC  } from '../../../modules/local/pave/somatic/main'

include { FileType                 } from '../utils_nfcore_oncoanalyser_pipeline/types'
include { groupByMeta              } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { joinMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { restoreMeta              } from '../utils_nfcore_oncoanalyser_pipeline/channel_helpers'
include { getInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getNormalDnaSample       } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSample        } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSampleName    } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { hasInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { hasNormalDna             } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { hasTumorDna              } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { selectCurrentOrExisting  } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow PAVE_ANNOTATION {
    take:
    // Sample data
    ch_inputs              // channel: [mandatory] [ meta ]
    ch_sage_dir_somatic    // channel: [mandatory] [ meta, sage_dir ]
    ch_sage_dir_germline   // channel: [mandatory] [ meta, sage_dir ]

    // Reference data
    genome_fasta           // channel: [mandatory] /path/to/genome_fasta
    genome_version         // channel: [mandatory] genome version
    genome_fai             // channel: [mandatory] /path/to/genome_fai
    pon_artefacts          // channel: [optional]  /path/to/pon_artefacts
    sage_pon               // channel: [mandatory] /path/to/sage_pon
    sage_blocklist_regions // channel: [mandatory] /path/to/sage_blocklist_regions
    sage_blocklist_sites   // channel: [mandatory] /path/to/sage_blocklist_sites
    clinvar_annotations    // channel: [mandatory] /path/to/clinvar_annotations
    segment_mappability    // channel: [mandatory] /path/to/segment_mappability
    driver_gene_panel      // channel: [mandatory] /path/to/driver_gene_panel
    ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/
    gnomad_resource        // channel: [mandatory] /path/to/gnomad_resource

    // Params
    sequencing_platform    // string:  [mandatory] sequencing platform

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
            def sage_vcf = sage_dir_selected ? sage_dir_selected.resolve("${getTumorDnaSampleName(meta)}.sage.germline.vcf.gz") : null
            def sage_tbi = sage_dir_selected ? sage_dir_selected.resolve("${getTumorDnaSampleName(meta)}.sage.germline.vcf.gz.tbi") : null

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

            def meta_pave = [
                key: meta.case_id,
                id: meta.case_id,
                sample_id: getTumorDnaSampleName(meta),
            ]

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
            def sage_vcf = sage_dir_selected ? sage_dir_selected.resolve("${getTumorDnaSampleName(meta)}.sage.somatic.vcf.gz") : null
            def sage_tbi = sage_dir_selected ? sage_dir_selected.resolve("${getTumorDnaSampleName(meta)}.sage.somatic.vcf.gz.tbi") : null

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

            def meta_pave = [
                key: meta.case_id,
                id: meta.case_id,
                sample_id: getTumorDnaSampleName(meta),
            ]

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
    // STEP: Outputs
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
