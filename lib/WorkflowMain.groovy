//
// This file holds several functions specific to the main.nf workflow in the nf-core/oncoanalyser pipeline
//

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            // TODO nf-core: Add Zenodo DOI for pipeline after first release
            //"* The pipeline\n" +
            //"  https://doi.org/10.5281/zenodo.XXXXXXX\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    //
    // Print help to screen if required
    //
    public static String help(workflow, params, log) {
        def command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome GRCh38 -profile docker"
        def help_string = ''
        help_string += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        help_string += NfcoreSchema.paramsHelp(workflow, params, command)
        help_string += '\n' + citation(workflow) + '\n'
        help_string += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return help_string
    }

    //
    // Print parameter summary log
    //
    public static String paramsSummaryLog(workflow, params, log) {
        def summary_log = ''
        summary_log += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        summary_log += NfcoreSchema.paramsSummaryLog(workflow, params)
        summary_log += '\n' + citation(workflow) + '\n'
        summary_log += NfcoreTemplate.dashedLine(params.monochrome_logs)
        log.info summary_log
    }

    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

        // Validate workflow parameters via the JSON schema
        if (params.validate_params) {
            NfcoreSchema.validateParameters(workflow, params, log)
        }

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)

        // Check that conda channels are set-up correctly
        if (params.enable_conda) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

        // Check input has been provided
        if (!params.input) {
            log.error "Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'"
            System.exit(1)
        }
    }
    //
    // Get attribute from genome config file e.g. fasta
    //
    public static Object getGenomeAttribute(params, attribute) {
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            if (params.genomes[ params.genome ].containsKey(attribute)) {
                return params.genomes[ params.genome ][ attribute ]
            }
        }
        return null
    }


    //
    // Set parameter defaults where required
    //
    public static void setParamsDefaults(params, log) {

        // Set defaults common to all run configuration

        def default_invalid = false

        if (params.containsKey('genome_version')) {
            params.ref_data_genome_version = params.genome_version.toString()
        } else if (Constants.GENOMES_VERSION_37.contains(params.genome)) {
            params.ref_data_genome_version = '37'
        } else if (Constants.GENOMES_VERSION_38.contains(params.genome)) {
            params.ref_data_genome_version = '38'
        } else {
            default_invalid = true
        }

        if (params.containsKey('genome_type')) {
            params.ref_data_genome_type = params.genome_type
        } else if (Constants.GENOMES_ALT.contains(params.genome)) {
            params.ref_data_genome_type = 'alt'
        } else if (Constants.GENOMES_DEFINED.contains(params.genome)) {
            params.ref_data_genome_type = 'no_alt'
        } else {
            default_invalid = true
        }

        if (!params.containsKey('ref_data_hmf_data_path')) {
            if (params.ref_data_genome_version == '37') {
                params.ref_data_hmf_data_path = Constants.HMF_DATA_37_PATH
            } else if (params.ref_data_genome_version == '38') {
                params.ref_data_hmf_data_path = Constants.HMF_DATA_38_PATH
            }
        }

        // Bad configuration, catch in validateParams
        if (default_invalid) {
            return
        }

        // Set defaults specific to run configuration without attempting to validate

        def run_mode
        if (params.containsKey('run_mode')) {
            run_mode = Utils.getRunMode(params.run_mode, log)
        } else {
            // Bad configuration, catch in validateParams
            return
        }

        if (run_mode == Constants.RunMode.RNA) {
            if (!params.containsKey('run_type')) {
                params.run_type = 'tumor_only'
            }
        }

        if (params.targeted === true) {

            if (!params.containsKey('run_type')) {
                params.run_type = 'tumor_only'
            }

            // Attempt to set default panel data path; make no assumption on valid 'panel' value
            if (!params.containsKey('ref_data_panel_data_path') && params.containsKey('panel')) {
                if (params.panel == 'hmf' && params.ref_data_genome_version == '38') {
                    params.ref_data_panel_data_path = Constants.HMF_PANEL_38_PATH
                } else if (params.panel == 'tso500' && params.ref_data_genome_version == '37') {
                    params.ref_data_panel_data_path = Constants.TSO500_PANEL_37_PATH
                } else if (params.panel == 'tso500' && params.ref_data_genome_version == '38') {
                    params.ref_data_panel_data_path = Constants.TSO500_PANEL_38_PATH
                }
            }

        }

        def stages = Processes.getRunStages(
            run_mode,
            params.processes_include,
            params.processes_exclude,
            params.processes_manual,
            params.targeted,
            log,
        )

        if (stages.virusinterpreter) {
            if (!params.containsKey('ref_data_virusbreakenddb_path')) {
                params.ref_data_virusbreakenddb_path = Constants.VIRUSBREAKENDDB_PATH
            }
        }

        if (stages.lilac) {
            if (params.ref_data_genome_version == '38' && params.ref_data_genome_type == 'alt' && !params.containsKey('ref_data_hla_slice_bed')) {
                params.ref_data_hla_slice_bed = Constants.HLA_SLICE_BED_GRCH38_ALT_PATH
            }
        }
    }

    //
    // Check and validate parameters
    //
    public static void validateParams(params, log) {

        // Common parameters

        if (!params.genome) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome must be set using the --genome CLI argument or in a configuration file.\n" +
                "  Currently, the available genome are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        } else if (!params.genomes.containsKey(params.genome)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        if (!Constants.GENOMES_SUPPORTED.contains(params.genome)) {
            if (!params.containsKey('force_genome') || !params.force_genome) {
                log.error "ERROR: currently only the GRCh37_hmf and GRCh38_hmf genomes are supported but got ${params.genome}" +
                    ", please adjust the --genome argument accordingly or override with --force_genome."
                System.exit(1)
            } else {
                log.warn "currently only the GRCh37_hmf and GRCh38_hmf genomes are supported but forcing to " +
                    "proceed with \"${params.genome}\""
            }
        }

        if (!params.ref_data_genome_version) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome version wasn't provided and genome '${params.genome}' is not defined in   \n" +
                "  genome version list.\n" +
                "  Currently, the list of genomes in the version list include:\n" +
                "  ${Constants.GENOMES_DEFINED.join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        if (!params.ref_data_genome_type) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome type wasn't provided and genome '${params.genome}' is not defined in      \n" +
                "  genome version list.\n" +
                "  Currently, the list of genomes in the version list include:\n" +
                "  ${Constants.GENOMES_DEINFED.join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        if (!params.ref_data_hmf_data_path) {
            log.error "ERROR: HMF data path wasn't provided"
            System.exit(1)
        }

        // NOTE(SW): this could be moved to the wgts.nf where we check that input files exist
        def null_check = [
           'ref_data_genome_fasta',
           'ref_data_genome_type',
           'ref_data_genome_version',
        ]
        null_check.each { k ->
            if (!params[k]) {
                log.error "ERROR: '${k}' cannot be set to null in any configuration and must be adjusted or removed to proceed."
                System.exit(1)
            }
        }

        // Run configuration specific parameters

        if (!params.run_mode) {
            def run_modes = Utils.getEnumNames(Constants.RunMode).join('\n    - ')
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Run mode must be set using the --run_mode CLI argument or in a configuration  \n" +
                "  file.\n" +
                "  Currently, the available run modes are:\n" +
                "    - ${run_modes}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        if (!params.run_type) {
            def run_types = Utils.getEnumNames(Constants.RunType).join('\n    - ')
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Run type must be set using the --run_type CLI argument or in a configuration  \n" +
                "  file.\n" +
                "  Currently, the available run types are:\n" +
                "    - ${run_types}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        def run_mode = Utils.getRunMode(params.run_mode, log)
        def run_type = Utils.getRunType(params.run_type, log)

        if (run_mode == Constants.RunMode.RNA) {
            if (run_type != Constants.RunType.TUMOR_ONLY) {
                log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  The WTS run mode does not support tumor/normal data, please adjust the CLI \n" +
                    "  --run_type option or corresponding configuration file.\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                System.exit(1)
            }

        }

        if (params.targeted === true) {

            if (run_type != Constants.RunType.TUMOR_ONLY) {
                log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  The panel run mode does not support tumor/normal data, please adjust the CLI \n" +
                    "  --run_type option or corresponding configuration file.\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                System.exit(1)
            }

            if (!params.containsKey('panel')) {

                def panels = Constants.PANELS_DEFINED.join('\n    - ')
                log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  A panel is required to be set using the --panel CLI argument or in a \n" +
                    "  configuration file.\n" +
                    "  Currently, the available panels are:\n" +
                    "    - ${panels}\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                System.exit(1)

            } else if (!Constants.PANELS_DEFINED.contains(params.panel)) {

                def panels = Constants.PANELS_DEFINED.join('\n    - ')
                log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  The ${params.panel} is not defined. Currently, the available panels are:\n" +
                    "    - ${panels}\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                System.exit(1)

            }

            if (params.panel == 'hmf' && params.ref_data_genome_version == '37') {
                log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  The Hartwig panel (hmf) is not available for the GRCh37 reference genome.\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                System.exit(1)
            }

            def rna_run = run_mode == Constants.RunMode.RNA || run_mode == Constants.RunMode.DNA_RNA
            if (rna_run && params.panel != 'tso500') {
                log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  Currently only the TSO500 panel (tso500) supports RNA analysis.\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                System.exit(1)
            }

        }

    }

    public static getRunConfig(params, log) {
        def run_mode = Utils.getRunMode(params.run_mode, log)
        def run_type = Utils.getRunType(params.run_type, log)
        def stages = Processes.getRunStages(
            run_mode,
            params.processes_include,
            params.processes_exclude,
            params.processes_manual,
            params.targeted,
            log,
        )

        return [
            mode: run_mode,
            type: run_type,
            stages: stages,
        ]
    }
}
