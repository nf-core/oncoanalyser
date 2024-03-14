//
// This file holds several functions specific to the main.nf workflow in the nf-core/oncoanalyser pipeline
//

class WorkflowMain {

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
        if (params.containsKey('mode') && params.mode !== null) {
            run_mode = Utils.getRunMode(params.mode, log)
        } else {
            // Bad configuration, catch in validateParams
            return
        }

        if (run_mode === Constants.RunMode.TARGETED) {

            // Attempt to set default panel data path; make no assumption on valid 'panel' value
            if (!params.containsKey('ref_data_panel_data_path') && params.containsKey('panel')) {
                if (params.panel == 'tso500' && params.ref_data_genome_version == '37') {
                    params.ref_data_panel_data_path = Constants.TSO500_PANEL_37_PATH
                } else if (params.panel == 'tso500' && params.ref_data_genome_version == '38') {
                    params.ref_data_panel_data_path = Constants.TSO500_PANEL_38_PATH
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
            if (!params.containsKey('ref_data_virusbreakenddb_path')) {
                params.ref_data_virusbreakenddb_path = Constants.VIRUSBREAKENDDB_PATH
            }
        }

        if (stages.lilac) {
            if (params.ref_data_genome_version == '38' && params.ref_data_genome_type == 'alt' && !params.containsKey('ref_data_hla_slice_bed')) {
                params.ref_data_hla_slice_bed = Constants.HLA_SLICE_BED_GRCH38_ALT_PATH
            }
        }

        if (stages.isofox && !params.containsKey('isofox_read_length')) {

            if (run_mode === Constants.RunMode.WGTS) {
                params.isofox_read_length = Constants.DEFAULT_ISOFOX_READ_LENGTH_WTS
            } else if (run_mode === Constants.RunMode.TARGETED) {
                params.isofox_read_length = Constants.DEFAULT_ISOFOX_READ_LENGTH_TARGETED
            }

        }

    }

    //
    // Check and validate parameters
    //
    public static void validateParams(params, log) {

        // Common parameters

        if (!params.genome) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
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
                "  genome type list.\n" +
                "  Currently, the list of genomes in the type list include:\n" +
                "  ${Constants.GENOMES_DEFINED.join(", ")}\n" +
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
