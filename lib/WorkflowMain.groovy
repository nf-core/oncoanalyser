//
// This file holds several functions specific to the main.nf workflow in the nf-core/oncoanalyser pipeline
//

class WorkflowMain {

    //
    // Set parameter defaults where required
    //
    public static void setParamsDefaults(params, log) {

        def default_invalid = false

        // Set defaults common to all run configuration

        if (params.genome_version !== null) {
            params.ref_data.genome_version = params.genome_version.toString()
        } else if (Constants.GENOMES_VERSION_37.contains(params.genome)) {
            params.ref_data.genome_version = '37'
        } else if (Constants.GENOMES_VERSION_38.contains(params.genome)) {
            params.ref_data.genome_version = '38'
        } else {
            default_invalid = true
        }

        if (params.genome_type !== null) {
            params.ref_data.genome_type = params.genome_type
        } else if (Constants.GENOMES_ALT.contains(params.genome)) {
            params.ref_data.genome_type = 'alt'
        } else if (Constants.GENOMES_DEFINED.contains(params.genome)) {
            params.ref_data.genome_type = 'no_alt'
        } else {
            default_invalid = true
        }

        if (params.hmf_data_path !== null) {
            params.ref_data.hmf_data_path = params.hmf_data_path
        } else if (params.ref_data.genome_version == '37') {
            params.ref_data.hmf_data_path = Constants.HMF_DATA_37_PATH
        } else if (params.ref_data.genome_version == '38') {
            params.ref_data.hmf_data_path = Constants.HMF_DATA_38_PATH
        } else {
            default_invalid = true
        }

        // Bad configuration, catch in validateParams
        if (default_invalid) {
            return
        }

        // Set defaults specific to run configuration without attempting to validate

        def run_mode
        if (params.mode !== null) {
            run_mode = Utils.getRunMode(params.mode, log)
        } else {
            // Bad configuration, catch in validateParams
            return
        }

        if (run_mode === Constants.RunMode.TARGETED) {

            // Attempt to set default panel data path; make no assumption on valid 'panel' value

            if (params.panel_data_path !== null) {
                params.ref_data.panel_data_path = params.panel_data_path
            } else if (params.panel !== null ) {
                if (params.panel == 'tso500' && params.genome_version == '37') {
                    params.ref_data.panel_data_path = Constants.TSO500_PANEL_37_PATH
                } else if (params.panel == 'tso500' && params.genome_version == '38') {
                    params.ref_data.panel_data_path = Constants.TSO500_PANEL_38_PATH
                }
            }
        }

        def stages = Processes.getRunStages(
            params.processes_include,
            params.processes_exclude,
            params.processes_manual,
            log,
        )

        if (stages.virusinterpreter && run_mode === Constants.RunMode.WGTS) {
            if (params.virusbreakenddb_path !== null) {
                params.ref_data.virusbreakenddb_path = params.virusbreakenddb_path
            } else {
                params.ref_data.virusbreakenddb_path = Constants.VIRUSBREAKENDDB_PATH
            }
        }

        if (stages.lilac) {
            if (params.hla_slice_bed !== null) {
                params.ref_data.hla_slice_bed = params.hla_slice_bed
            } else if (params.genome_version == '38' && params.genome_type == 'alt') {
                params.ref_data.hla_slice_bed = Constants.HLA_SLICE_BED_GRCH38_ALT_PATH
            }
        }

    }

    //
    // Check and validate parameters
    //
    public static void validateParams(params, log) {

        // Common parameters

        if (!params.ref_data.genome) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome must be set using the --genome CLI argument or in a configuration file.\n" +
                "  Currently, the available genome are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        } else if (!params.genomes.containsKey(params.ref_data.genome)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.ref_data.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        if (!Constants.GENOMES_SUPPORTED.contains(params.ref_data.genome)) {
            if (!params.ref_data.force_genome) {
                log.error "ERROR: currently only the GRCh37_hmf and GRCh38_hmf genomes are supported but got ${params.ref_data.genome}" +
                    ", please adjust the --genome argument accordingly or override with --force_genome."
                System.exit(1)
            } else {
                log.warn "currently only the GRCh37_hmf and GRCh38_hmf genomes are supported but forcing to " +
                    "proceed with \"${params.ref_data.genome}\""
            }
        }

        if (!params.ref_data.genome_version) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome version wasn't provided and genome '${params.ref_data.genome}' is not defined in   \n" +
                "  genome version list.\n" +
                "  Currently, the list of genomes in the version list include:\n" +
                "  ${Constants.GENOMES_DEFINED.join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        if (!params.ref_data.genome_type) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome type wasn't provided and genome '${params.ref_data.genome}' is not defined in      \n" +
                "  genome type list.\n" +
                "  Currently, the list of genomes in the type list include:\n" +
                "  ${Constants.GENOMES_DEFINED.join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        if (!params.ref_data.hmf_data_path) {
            log.error "ERROR: HMF data path wasn't provided"
            System.exit(1)
        }

        // Run configuration specific parameters

        if (!params.mode) {
            def run_modes = Utils.getEnumNames(Constants.RunMode).join('\n    - ')
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Run mode must be set using the --mode CLI argument or in a configuration  \n" +
                "  file.\n" +
                "  Currently, the available run modes are:\n" +
                "    - ${run_modes}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        def run_mode = Utils.getRunMode(params.mode, log)

        if (run_mode === Constants.RunMode.TARGETED) {

            if (!params.containsKey('panel')) {

                def panels = Constants.PANELS_DEFINED.join('\n    - ')
                log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  A panel is required to be set using the --panel CLI argument or in a \n" +
                    "  configuration file.\n" +
                    "  Currently, the available panels are:\n" +
                    "    - ${panels}\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                System.exit(1)

            } else if (!Constants.PANELS_DEFINED.contains(params.panel)) {

                def panels = Constants.PANELS_DEFINED.join('\n    - ')
                log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  The ${params.panel} is not defined. Currently, the available panels are:\n" +
                    "    - ${panels}\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                System.exit(1)

            }

        }

    }

    public static getRunConfig(params, inputs, log) {

        def run_mode = Utils.getRunMode(params.mode, log)

        def stages = Processes.getRunStages(
            params.processes_include,
            params.processes_exclude,
            params.processes_manual,
            log,
        )

        return [
            mode: run_mode,
            panel: run_mode === Constants.RunMode.TARGETED ? params.panel : null,
            stages: stages,
            has_dna: inputs.any { it.containsKey([Constants.SampleType.TUMOR, Constants.SequenceType.DNA]) },
            has_rna: inputs.any { it.containsKey([Constants.SampleType.TUMOR, Constants.SequenceType.RNA]) },
        ]
    }
}
