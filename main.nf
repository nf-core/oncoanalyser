#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/oncoanalyser
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/oncoanalyser
    Website: https://nf-co.re/oncoanalyser
    Slack  : https://nfcore.slack.com/channels/oncoanalyser
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.types = true

include { PIPELINE_INITIALISATION  } from './subworkflows/local/utils_nfcore_oncoanalyser_pipeline'
include { PIPELINE_COMPLETION      } from './subworkflows/local/utils_nfcore_oncoanalyser_pipeline'

include { getGenomeAttribute       } from './subworkflows/local/utils_nfcore_oncoanalyser_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SET GENOME VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {
    ref_data_genome_fasta         = getGenomeAttribute('fasta')
    ref_data_genome_fai           = getGenomeAttribute('fai')
    ref_data_genome_dict          = getGenomeAttribute('dict')
    ref_data_genome_img           = getGenomeAttribute('img')
    ref_data_genome_bwamem2_index = getGenomeAttribute('bwamem2_index')
    ref_data_genome_gridss_index  = getGenomeAttribute('gridss_index')
    ref_data_genome_star_index    = getGenomeAttribute('star_index')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PANEL_RESOURCE_CREATION  } from './workflows/panel_resource_creation'
include { PREPARE_REFERENCE        } from './workflows/prepare_reference'
include { PURITY_ESTIMATE          } from './workflows/purity_estimate'
include { TARGETED                 } from './workflows/targeted'
include { WGTS                     } from './workflows/wgts'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { parseInput             } from './subworkflows/local/utils_nfcore_oncoanalyser_pipeline/parse_inputs'
include { RunMode                } from './subworkflows/local/utils_nfcore_oncoanalyser_pipeline/types'
include { createStubPlaceholders } from './subworkflows/local/utils_nfcore_oncoanalyser_pipeline/utils'
include { getRunMode             } from './subworkflows/local/utils_nfcore_oncoanalyser_pipeline/utils'
include { getRunConfig           } from './subworkflows/local/utils_nfcore_oncoanalyser_pipeline/validate_params'
include { setParamsDefaults      } from './subworkflows/local/utils_nfcore_oncoanalyser_pipeline/validate_params'
include { validateInput          } from './subworkflows/local/utils_nfcore_oncoanalyser_pipeline/validate_params'
include { validateParams         } from './subworkflows/local/utils_nfcore_oncoanalyser_pipeline/validate_params'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//

workflow NFCORE_ONCOANALYSER {
    main:
    // Get run mode
    def run_mode = getRunMode(params.mode, log)

    // Results channel for eventual publishing
    // channel: [filepath, file]
    ch_results = channel.empty()

    // Run selected workflow
    // NOTE(SW): prepare reference is checked early as params.input is not required
    if (run_mode == RunMode.PREPARE_REFERENCE)  {
        PREPARE_REFERENCE(params)
        ch_results = PREPARE_REFERENCE.out.results
    } else {
        // Parse and validate inputs
        inputs = parseInput(params.input, workflow.stubRun, log)
        run_config = getRunConfig(params, inputs, log)
        validateInput(inputs, run_config, params, log)

        // Run requested workflow
        if (run_mode == RunMode.WGTS) {
            WGTS(inputs, run_config, params)
            ch_results = WGTS.out.results
        } else if (run_mode == RunMode.TARGETED) {
            TARGETED(inputs, run_config, params)
            ch_results = TARGETED.out.results
        } else if (run_mode == RunMode.PURITY_ESTIMATE) {
            PURITY_ESTIMATE(inputs, run_config, params)
            ch_results = PURITY_ESTIMATE.out.results
        } else if (run_mode == RunMode.PANEL_RESOURCE_CREATION) {
            PANEL_RESOURCE_CREATION(inputs, run_config, params)
            ch_results = PANEL_RESOURCE_CREATION.out.results
        } else {
            log.error("received bad run mode: ${run_mode}")
            exit(1)
        }
    }

    emit:
    results = ch_results // channel: [ filepath, file ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    main:
    //
    // STEP: Set defaults and apply extended, custom validation
    //
    setParamsDefaults(params, log)
    validateParams(params, log)

    //
    // STEP: Create placeholders for stub runs if requested
    //
    if (workflow.stubRun && params.create_stub_placeholders) {
        createStubPlaceholders(params)
    }

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_ONCOANALYSER()

    ch_results = NFCORE_ONCOANALYSER.out.results
    // NOTE(SW): extract the MultiQC report from results instead of re-reading the
    // single-consumer 'multiqc_report' topic (already read by multiqc_reporting)
    ch_multiqc_report = ch_results
        .filter { filepath, file -> file.name == 'multiqc_report.html' }
        .map { filepath, file -> file }

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        ch_multiqc_report,
    )

    publish:
    results = ch_results
}

output {
    results {
        path { filepath, file -> file >> filepath }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
