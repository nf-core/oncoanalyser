/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.types = true

include { PREPARE_OUTPUTS           } from '../subworkflows/local/prepare_outputs'
include { PREPARE_REFERENCE as STAGE_REFERENCE } from '../subworkflows/local/prepare_reference'

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { getPrepConfigFromCli } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/validate_params'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPARE_REFERENCE {
    take:
    params

    main:
    // Stage in reference data as requested
    def prep_config = getPrepConfigFromCli(params, log)
    STAGE_REFERENCE(
        prep_config,
        [:],
        params,
    )

    //
    // TASK: Aggregate software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'oncoanalyser_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true,
        )

    //
    // SUBWORKFLOW: Prepare outputs for publishing
    //
    PREPARE_OUTPUTS(
        channel.empty(),  // amber
        channel.empty(),  // bamtools_tumor
        channel.empty(),  // bamtools_normal
        channel.empty(),  // chord
        channel.empty(),  // cider
        channel.empty(),  // cobalt
        channel.empty(),  // cuppa
        channel.empty(),  // esvee
        channel.empty(),  // align_rna_tumor
        channel.empty(),  // isofox
        channel.empty(),  // lilac
        channel.empty(),  // linx_germline
        channel.empty(),  // linx_somatic
        channel.empty(),  // linx_somatic_visualiser
        channel.empty(),  // linxreport_html
        channel.empty(),  // multiqc
        channel.empty(),  // neo_annotated_fusions
        channel.empty(),  // neo_finder
        channel.empty(),  // neo_scorer_dir
        channel.empty(),  // orange_json
        channel.empty(),  // orange_pdf
        channel.empty(),  // pave_germline
        channel.empty(),  // pave_somatic
        channel.empty(),  // peach
        channel.empty(),  // purple
        channel.empty(),  // qsee
        channel.empty(),  // redux_tumor
        channel.empty(),  // redux_normal
        channel.empty(),  // redux_donor
        channel.empty(),  // sage_append_somatic
        channel.empty(),  // sage_append_germline
        channel.empty(),  // sage_germline
        channel.empty(),  // sage_somatic
        channel.empty(),  // sage_somatic_visualiser
        channel.empty(),  // sigs
        channel.empty(),  // teal_normal_bam
        channel.empty(),  // teal_tumor_bam
        channel.empty(),  // teal_tsvs
        channel.empty(),  // virusbreakend_tsv
        channel.empty(),  // virusbreakend_vcf
        channel.empty(),  // virusinterpreter
        channel.topic('write_reference_data'),
        channel.topic('command_files'),
        channel.empty(),  // wisp
        channel.empty(),  // cobalt_normalisation_tsv
        channel.empty(),  // isofox_normalisation_csv
        channel.empty(),  // pave_pon_panel_creation_artefacts
    )

    emit:
    results = PREPARE_OUTPUTS.out.results
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
