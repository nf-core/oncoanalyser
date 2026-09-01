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

include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_oncoanalyser_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_oncoanalyser_pipeline'

include { getGenomeAttribute } from './subworkflows/local/utils_nfcore_oncoanalyser_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SET GENOME VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.ref_data_genome_fasta         = getGenomeAttribute('fasta')
params.ref_data_genome_fai           = getGenomeAttribute('fai')
params.ref_data_genome_dict          = getGenomeAttribute('dict')
params.ref_data_genome_img           = getGenomeAttribute('img')
params.ref_data_genome_bwamem2_index = getGenomeAttribute('bwamem2_index')
params.ref_data_genome_gridss_index  = getGenomeAttribute('gridss_index')
params.ref_data_genome_star_index    = getGenomeAttribute('star_index')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PANEL_RESOURCE_CREATION } from './workflows/panel_resource_creation'
include { PREPARE_REFERENCE       } from './workflows/prepare_reference'
include { PURITY_ESTIMATE         } from './workflows/purity_estimate'
include { TARGETED                } from './workflows/targeted'
include { WGTS                    } from './workflows/wgts'

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
    def run_mode = Utils.getRunMode(params.mode, log)

    // Results channel for eventual publishing
    // channel: [filepath, file]
    ch_results = channel.empty()

    // Run selected workflow
    // NOTE(SW): prepare reference is checked early as params.input is not required
    if (run_mode == Constants.RunMode.PREPARE_REFERENCE)  {
        PREPARE_REFERENCE(params)
        ch_results = ch_results.mix(PREPARE_REFERENCE.out.results)
    } else {
        // Parse and validate inputs
        inputs = Utils.parseInput(params.input, workflow.stubRun, log)
        run_config = WorkflowMain.getRunConfig(params, inputs, log)
        Utils.validateInput(inputs, run_config, params, log)

        // Run requested workflow
        if (run_mode == Constants.RunMode.WGTS) {
            WGTS(inputs, run_config, params)
            ch_results = ch_results.mix(WGTS.out.results)
        } else if (run_mode == Constants.RunMode.TARGETED) {
            TARGETED(inputs, run_config, params)
            ch_results = ch_results.mix(TARGETED.out.results)
        } else if (run_mode == Constants.RunMode.PURITY_ESTIMATE) {
            PURITY_ESTIMATE(inputs, run_config, params)
            ch_results = ch_results.mix(PURITY_ESTIMATE.out.results)
        } else if (run_mode == Constants.RunMode.PANEL_RESOURCE_CREATION) {
            PANEL_RESOURCE_CREATION(inputs, run_config, params)
            ch_results = ch_results.mix(PANEL_RESOURCE_CREATION.out.results)
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
    WorkflowMain.setParamsDefaults(params, log)
    WorkflowMain.validateParams(params, log)

    //
    // STEP: Create placeholders for stub runs if requested
    //
    if (workflow.stubRun && params.create_stub_placeholders) {
        Utils.createStubPlaceholders(params)
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

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        channel.topic('multiqc_report'),
    )

    publish:
    results = NFCORE_ONCOANALYSER.out.results
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
