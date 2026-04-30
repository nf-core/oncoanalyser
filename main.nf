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
    SET DEFAULT VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.ref_data_genome_fasta         = getGenomeAttribute('fasta')
params.ref_data_genome_fai           = getGenomeAttribute('fai')
params.ref_data_genome_dict          = getGenomeAttribute('dict')
params.ref_data_genome_img           = getGenomeAttribute('img')
params.ref_data_genome_bwamem2_index = getGenomeAttribute('bwamem2_index')
params.ref_data_genome_gridss_index  = getGenomeAttribute('gridss_index')
params.ref_data_genome_star_index    = getGenomeAttribute('star_index')

pipeline.Params.parse(params)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE PLACEHOLDER FILES FOR STUB RUNS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// NOTE(SW): required prior to workflow import

if (workflow.stubRun && params.create_stub_placeholders) {
    pipeline.Params.createStubPlaceholders(params)
}

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

    def pipeline_mode = pipeline.PipelineMode.fromString(params.mode)

    // Run selected workflow
    // NOTE(SW): prepare reference is checked early as params.input is not required
    if (pipeline_mode === pipeline.PipelineMode.PREPARE_REFERENCE)  {
        PREPARE_REFERENCE()
    } else {

        def inputs = samplesheet.SampleSheet.parse(params.input, pipeline_mode)

        def stages = pipeline.RunStage.getValidatedRunStages(
           params.processes_include,
           params.processes_exclude,
           params.processes_manual,
           log,
        )

        // Run requested workflow
        if (pipeline_mode === pipeline.PipelineMode.WGTS) {
            WGTS(inputs, stages)
        } else if (pipeline_mode === pipeline.PipelineMode.TARGETED) {
            TARGETED(inputs, stages)
        } else if (pipeline_mode === pipeline.PipelineMode.PURITY_ESTIMATE) {
            PURITY_ESTIMATE(inputs, stages)
        } else if (pipeline_mode === pipeline.PipelineMode.PANEL_RESOURCE_CREATION) {
            PANEL_RESOURCE_CREATION(inputs, stages)
        } else {
            log.error("received bad pipeline mode: ${pipeline_mode}")
            Nextflow.exit(1)
        }
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

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
        params.hook_url,
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
