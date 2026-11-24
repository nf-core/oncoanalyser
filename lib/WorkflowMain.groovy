//
// This file holds several functions specific to the main.nf workflow in the nf-core/oncoanalyser pipeline
//

import nextflow.Nextflow

import Utils

class WorkflowMain {

    //
    // Set parameter defaults where required
    //
    public static void setParamsDefaults(params, log) {

        def default_invalid = false

        // Set defaults common to all run configuration
        if (! params.containsKey('sequencing_platform')) {
            params.sequencing_platform = 'illumina'
        }

        if (! params.containsKey('genome_version')) {
            if (Constants.GENOMES_VERSION_37.contains(params.genome)) {
                params.genome_version = '37'
            } else if (Constants.GENOMES_VERSION_38.contains(params.genome)) {
                params.genome_version = '38'
            } else {
                default_invalid = true
            }
        }

        if (! params.containsKey('genome_type')) {
            if (Constants.GENOMES_ALT.contains(params.genome)) {
                params.genome_type = 'alt'
            } else if (Constants.GENOMES_DEFINED.contains(params.genome)) {
                params.genome_type = 'no_alt'
            } else {
                default_invalid = true
            }
        }

        if (! params.containsKey('ref_data_hmf_data_path')) {
            if (params.genome_version.toString() == '37') {
                params.ref_data_hmf_data_path = "${params.ref_data_base}/${Constants.HMF_DATA_37_PATH}"
            } else if (params.genome_version.toString() == '38') {
                params.ref_data_hmf_data_path = "${params.ref_data_base}/${Constants.HMF_DATA_38_PATH}"
            } else {
                default_invalid = true
            }
        }

        // Bad configuration, catch in validateParams
        if (default_invalid) {
            return
        }

        // Set defaults specific to run configuration without attempting to validate

        def run_mode
        if (params.mode != null) {
            run_mode = Utils.getRunMode(params.mode, log)
        } else {
            // Bad configuration, catch in validateParams
            return
        }

        // Attempt to set default panel data path; make no assumption on valid 'panel' value
        if (run_mode == Constants.RunMode.TARGETED || run_mode == Constants.RunMode.PREPARE_REFERENCE) {

            if (params.containsKey('panel')) {

                if (params.panel.toLowerCase() == 'tso500') {
                    if (params.genome_version.toString() == '37') {
                        params.ref_data_panel_data_path = "${params.ref_data_base}/${Constants.TSO500_PANEL_37_PATH}"
                    }
                }

            }

        }

        def stages = Processes.getRunStages(
            params.processes_include,
            params.processes_exclude,
            params.processes_manual,
            log,
        )

        //
        // Resolve UMI type and set UMI parameters
        //
        def umi_type
        if (params.containsKey('umi_type') && params.umi_type) {
            umi_type = Utils.getEnumFromString(params.umi_type, Constants.UmiType)
        } else if (params.containsKey('panel') && Constants.PANELS_DEFINED.contains(params.panel.toLowerCase())) {
            if (params.panel.toLowerCase() == 'tso500') {
                umi_type = Constants.UmiType.TSO500
            }
        }

        if (umi_type == Constants.UmiType.KAPA) {
              params.fastp_umi_enabled = true
              params.fastp_umi_location = 'per_read'
              params.fastp_umi_length = 6
              params.fastp_umi_skip = 0
              params.redux_umi_enabled = true
              params.redux_umi_duplex_delim = '_'
        } else if (umi_type == Constants.UmiType.MSK) {
              params.fastq_tools_umi_enabled = true
              params.fastq_tools_umi_delim = '+'
              params.redux_umi_enabled = true
              params.redux_umi_duplex_delim = ''
        } else if (umi_type == Constants.UmiType.TSO500) {
              params.redux_umi_enabled = true
              params.redux_umi_duplex_delim = '+'
        } else if (umi_type == Constants.UmiType.TWIST) {
              params.fastp_umi_enabled = true
              params.fastp_umi_location = 'per_read'
              params.fastp_umi_length = 7
              params.fastp_umi_skip = 0
              params.redux_umi_enabled = true
              params.redux_umi_duplex_delim = '_'
        }

        // Final point to set any default to avoid access to undefined parameters during nf-validation
        if (! params.containsKey('panel')) params.panel = null
        if (! params.containsKey('ref_data_genome_alt')) params.ref_data_genome_alt = null
        if (! params.containsKey('ref_data_genome_gtf')) params.ref_data_genome_gtf = null
        if (! params.containsKey('ref_data_panel_data_path')) params.ref_data_panel_data_path = null

        // Additionally set selected parameters with false-ish truthy values to avoid passing null values as inputs
        if (! params.containsKey('fastp_umi_enabled')) params.fastp_umi_enabled = false
        if (! params.containsKey('fastp_umi_length')) params.fastp_umi_length = 0
        if (! params.containsKey('fastp_umi_location')) params.fastp_umi_location = ''
        if (! params.containsKey('fastp_umi_skip')) params.fastp_umi_skip = -1
        if (! params.containsKey('fastq_tools_umi_enabled')) params.fastq_tools_umi_enabled = false
        if (! params.containsKey('fastq_tools_umi_delim')) params.fastq_tools_umi_delim = ''
        if (! params.containsKey('redux_umi_enabled')) params.redux_umi_enabled = false
        if (! params.containsKey('redux_umi_duplex_delim')) params.redux_umi_duplex_delim = ''

    }

    //
    // Check and validate parameters
    //
    public static void validateParams(params, log) {

        // Common parameters

        if (! params.genome) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome must be set using the --genome CLI argument or in a configuration file.\n" +
                "  Currently, the available genomes are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        } else if (! params.genomes.containsKey(params.genome)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genomes are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        if (! Constants.GENOMES_SUPPORTED.contains(params.genome)) {
            if (! params.force_genome) {
                log.error "currently only the GRCh37_hmf and GRCh38_hmf genomes are supported but got ${params.genome}" +
                    ", please adjust the --genome argument accordingly or override with --force_genome."
                Nextflow.exit(1)
            } else {
                log.warn "currently only the GRCh37_hmf and GRCh38_hmf genomes are supported but forcing to " +
                    "proceed with \"${params.genome}\""
            }
        }

        if (! params.genome_version) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome version wasn't provided and genome '${params.genome}' is not defined in   \n" +
                "  genome version list.\n" +
                "  Currently, the list of genomes in the version list includes:\n" +
                "  ${Constants.GENOMES_DEFINED.join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        if (! params.genome_type) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome type wasn't provided and genome '${params.genome}' is not defined in      \n" +
                "  genome type list.\n" +
                "  Currently, the list of genomes in the type list include:\n" +
                "  ${Constants.GENOMES_DEFINED.join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        if (! params.ref_data_hmf_data_path) {
            log.error "HMF data path wasn't provided"
            Nextflow.exit(1)
        }

        if (! params.hmf_data_paths.containsKey(params.genome_version.toString())) {
            def hmf_data_versions = params.hmf_data_paths.keySet().join('\n    - ')
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Did not find path definitions for the provided 'hmf_data_paths' and\n" +
                "  genome version ${params.genome_version}. Please check your configuration.\n" +
                "  Found the following genome version definitions for 'hmf_data_paths':\n" +
                "    - ${hmf_data_versions}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        // Run configuration specific parameters

        if (! params.mode) {
            def run_modes = Utils.getEnumNames(Constants.RunMode).join('\n    - ')
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Run mode must be set using the --mode CLI argument or in a configuration  \n" +
                "  file.\n" +
                "  Currently, the available run modes are:\n" +
                "    - ${run_modes}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        def run_mode = Utils.getRunMode(params.mode, log)

        if (run_mode == Constants.RunMode.PREPARE_REFERENCE && params.ref_data_types == null) {

            def ref_data_types = Utils.getEnumNames(Constants.RefDataType).join('\n    - ')

            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  CLI argument --ref_data_types is required for mode prepare_reference.\n" +
                "  Please specify one or more of the below valid values (separated by commas)\n" +
                "    - ${ref_data_types}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        if (run_mode == Constants.RunMode.WGTS) {

            // We allow fastq-tools UMI processing in WGTS but this requires that the user manually set the 'known_umis' file in hmf_data, enforce here
            def has_known_umis = params.hmf_data_paths[params.genome_version.toString()].containsKey('known_umis')
            if (params.fastq_tools_umi_enabled && ! has_known_umis) {
                log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  The fastq-tools process is enabled but the required 'known_umis' data path was\n" +
                    "  not configured in the respective 'hmf_data_paths' for genome version ${params.genome_version.toString()}\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                Nextflow.exit(1)
            }

        }

        if (run_mode == Constants.RunMode.TARGETED) {

            if (! params.containsKey('panel') || params.panel == null) {

                def panels = Constants.PANELS_DEFINED.join('\n    - ')
                log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  A panel is required to be set using the --panel CLI argument or in a\n" +
                    "  configuration file when running in targeted mode or panel resource creation mode.\n" +
                    "  Currently, the available built-in panels are:\n" +
                    "    - ${panels}\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                Nextflow.exit(1)

            } else if (! Constants.PANELS_DEFINED.contains(params.panel.toLowerCase())) {

                if (params.containsKey('force_panel') && params.force_panel) {
                    log.warn "provided panel ${params.panel.toLowerCase()} does not have built-in support but forcing to proceed"
                } else {
                    def panels = Constants.PANELS_DEFINED.join('\n    - ')
                    log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                        "  The ${params.panel.toLowerCase()} panel does not have built-in support. Currently, the\n" +
                        "  available supported panels are:\n" +
                        "    - ${panels}\n\n" +
                        "  Please adjust the --panel argument or override with --force_panel.\n" +
                        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                    Nextflow.exit(1)
                }

            }

            // Require the panel to have defined

            if (! params.panel_data_paths.containsKey(params.panel.toLowerCase())) {
                def panels = params.panel_data_paths.keySet().join('\n    - ')
                log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  Did not find data path definitions for the provided ${params.panel.toLowerCase()} panel.\n" +
                    "  Please check your configuration. Found the following panel definitions:\n" +
                    "    - ${panels}\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                Nextflow.exit(1)
            }

            def panel_data_paths_versions = params.panel_data_paths[params.panel.toLowerCase()]
            if (! panel_data_paths_versions.containsKey(params.genome_version.toString())) {
                def panel_versions = panel_data_paths_versions.keySet().join('\n    - ')
                log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  Did not find path definitions for the provided ${params.panel.toLowerCase()} panel and\n" +
                    "  genome version ${params.genome_version}. Please check your configuration.\n" +
                    "  Found the following genome version panel definitions for ${params.panel.toLowerCase()}:\n" +
                    "    - ${panel_versions}\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                Nextflow.exit(1)
            }

            // Perform check for required and optional fields

            def panel_data_paths = panel_data_paths_versions[params.genome_version.toString()]

            def required_entries = [
                'driver_gene_panel',
                'pon_artefacts',
                'target_regions_bed',
                'target_regions_normalisation',
            ]

            def optional_entries = [
                'isofox_counts',
                'isofox_gc_ratios',
                'isofox_tpm_norm',
                'known_umis',
                'msi_model_error_rates',
            ]

            def required_entries_missing = required_entries.findAll { n -> ! panel_data_paths.containsKey(n) || ! panel_data_paths[n] }
            if (required_entries_missing) {
                def required_entries_missing_str = required_entries_missing.join('\n    - ')
                log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  The following panel data path entries are required but were not found:\n" +
                    "    - ${required_entries_missing_str}\n\n" +
                    "  Please review configuration for the ${params.panel.toLowerCase()} (${params.genome_version}) panel\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                Nextflow.exit(1)
            }

            // Require optional entries to also always be present, but either set to a file or an empty list []; i.e. do not allow empty string ''
            def optional_entries_missing = optional_entries.findAll { n -> ! panel_data_paths.containsKey(n) }
            if (optional_entries_missing) {
                def optional_entries_missing_str = optional_entries_missing.join('\n    - ')
                log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  The following panel data path entries are required but were not found:\n" +
                    "    - ${optional_entries_missing_str}\n\n" +
                    "  Please review configuration for the ${params.panel.toLowerCase()} (${params.genome_version}) panel\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                Nextflow.exit(1)
            }

            def optional_entries_invalid = optional_entries.findAll { n -> ! panel_data_paths[n] && panel_data_paths[n] != []  }
            if (optional_entries_invalid) {
                def optional_entries_invalid_str = optional_entries_invalid.join('\n    - ')
                log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  The following panel data path optional entries should only be set to [] if\n" +
                    "  not applicable:\n" +
                    "    - ${optional_entries_invalid_str}\n\n" +
                    "  Please review configuration for the ${params.panel.toLowerCase()} (${params.genome_version}) panel\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                Nextflow.exit(1)
            }

            // Require known_umis in panel data when fastq-tools UMI processing is enabled
            if (params.fastq_tools_umi_enabled && ! panel_data_paths['known_umis']) {
                log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  The fastq-tools process is enabled but the required 'known_umis' data path was\n" +
                    "  not configured in the panel data paths for panel ${params.panel.toLowerCase()}\n" +
                    "  (${params.genome_version}).\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                Nextflow.exit(1)
            }

        }

        if (run_mode == Constants.RunMode.PURITY_ESTIMATE) {

            def purity_estimate_modes = [Constants.RunMode.WGTS, Constants.RunMode.TARGETED]

            def purity_mode_enum = ! params.purity_estimate_mode ? null : Utils.getEnumFromString(params.purity_estimate_mode, Constants.RunMode)

            if (! purity_mode_enum || ! purity_estimate_modes.contains(purity_mode_enum)) {

                def purity_estimate_modes_str = purity_estimate_modes
                    .collect { e -> e.name().toLowerCase() }
                    .join('\n    - ')

                log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  A valid purity estimate run mode must be set using the --purity_estimate_mode\n" +
                    "  CLI argument or in a configuration file.\n" +
                    "  Currently, the available run modes are:\n" +
                    "    - ${purity_estimate_modes_str}\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                Nextflow.exit(1)
            }
        }

        if (params.ref_data_genome_alt != null) {
            if (params.genome_type != 'alt') {
                log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  Using a reference genome without ALT contigs but found an .alt file\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                Nextflow.exit(1)
            }

            def ref_data_genome_alt_fn = nextflow.Nextflow.file(params.ref_data_genome_alt).name
            def ref_data_genome_fasta_fn = nextflow.Nextflow.file(params.ref_data_genome_fasta).name
            if (ref_data_genome_alt_fn != "${ref_data_genome_fasta_fn}.alt") {
                log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  Found .alt file with filename of ${ref_data_genome_alt_fn} but it is required to match\n" +
                    "  reference genome FASTA filename stem: ${ref_data_genome_fasta_fn}.alt\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                Nextflow.exit(1)
            }

        }

        // UMI parameters
        if (params.fastp_umi_enabled && params.fastq_tools_umi_enabled) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                  "  UMI processing with either fastp or fastq-tools but not both can be enabled\n" +
                  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        if ((params.fastp_umi_enabled || params.fastq_tools_umi_enabled) && ! params.redux_umi_enabled) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                  "  When FASTQ UMI processing is enabled (via fastp_umi_enabled or fastq_tools_umi_enabled),\n" +
                  "  REDUX UMI processing must also be enabled with redux_umi_enabled\n" +
                  "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        def fastp_umi_args_set_any = params.fastp_umi_location || params.fastp_umi_length || params.fastp_umi_skip >= 0
        if (fastp_umi_args_set_any && ! params.fastp_umi_enabled) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Detected use of fastp UMI parameters but fastp UMI processing has not been enabled.\n" +
                "  Please review your configuration and set the fastp_umi_enabled flag or otherwise " +
                "  adjust accordingly.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        def fastp_umi_args_set_all = params.fastp_umi_location && params.fastp_umi_length && params.fastp_umi_skip >= 0
        if (params.fastp_umi_enabled && ! fastp_umi_args_set_all) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Refusing to run fastp UMI processing without having all UMI params configured.\n" +
                "  Please review your configuration and appropriately set all fastp_umi_* parameters.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        if (params.fastq_tools_umi_delim && ! params.fastq_tools_umi_enabled) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Detected use of fastq-tools UMI parameter fastq_tools_umi_delim but fastq-tools UMI\n" +
                "  processing has not been enabled. Please review your configuration and set the\n" +
                "  fastq_tools_umi_enabled or otherwise adjust accordingly.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        if (params.fastq_tools_umi_enabled && ! params.fastq_tools_umi_delim) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Refusing to run fastq-tools UMI processing without fastq_tools_umi_delim configured.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        if (params.redux_umi_duplex_delim && params.redux_umi_enabled == false) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Detected use of REDUX UMI parameters but REDUX UMI processing has not been\n" +
                "  enabled. Please review your configuration and set the redux_umi_enabled flag or\n" +
                "  otherwise adjust accordingly.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
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
            stages: stages,
            has_dna: inputs.any { Utils.hasTumorDna(it) },
            has_rna: inputs.any { Utils.hasTumorRna(it) },
            has_rna_fastq: inputs.any { Utils.hasTumorRnaFastq(it) },
            has_dna_fastq: inputs.any { Utils.hasTumorDnaFastq(it) || Utils.hasNormalDnaFastq(it) },
        ]
    }

    public static getPrepConfigFromSamplesheet(run_config) {
        return [
            prepare_ref_data_only: false,

            require_fasta: true,
            require_fai: true,
            require_dict: true,
            require_img: true,

            require_bwamem2_index: run_config.has_dna_fastq && run_config.stages.alignment,
            require_star_index: run_config.has_rna_fastq && run_config.stages.alignment,

            require_gridss_index: run_config.has_dna && run_config.mode == Constants.RunMode.WGTS && run_config.stages.virusinterpreter,
            require_hmftools_data: true,
            require_panel_data: run_config.mode == Constants.RunMode.TARGETED,
        ]
    }

    public static getPrepConfigFromCli(params, log) {
        def ref_data_types = params.ref_data_types
            .tokenize(',')
            .collect {
                def ref_data_type_enum = Utils.getEnumFromString(it, Constants.RefDataType)

                if (! ref_data_type_enum) {
                    def ref_data_type_str = Utils.getEnumNames(Constants.RefDataType).join('\n  - ')
                    log.error "received invalid ref data type: '${it}'. Valid options are:\n  - ${ref_data_type_str}"
                    Nextflow.exit(1)
                }

                return ref_data_type_enum
            }

        if (
            ref_data_types.contains(Constants.RefDataType.WGS) ||
            ref_data_types.contains(Constants.RefDataType.WTS) ||
            ref_data_types.contains(Constants.RefDataType.TARGETED)
        ) {
            ref_data_types += [
                Constants.RefDataType.FASTA,
                Constants.RefDataType.FAI,
                Constants.RefDataType.DICT,
                Constants.RefDataType.IMG,
                Constants.RefDataType.HMFTOOLS
            ]
        }

        if (ref_data_types.contains(Constants.RefDataType.WGS)) {
            ref_data_types += [Constants.RefDataType.GRIDSS_INDEX]
        }

        if (ref_data_types.contains(Constants.RefDataType.TARGETED)) {
            ref_data_types += [Constants.RefDataType.PANEL]
        }

        def require_fasta = ref_data_types.contains(Constants.RefDataType.FASTA)
        def require_fai = ref_data_types.contains(Constants.RefDataType.FAI)
        def require_dict = ref_data_types.contains(Constants.RefDataType.DICT)
        def require_img = ref_data_types.contains(Constants.RefDataType.IMG)

        def require_bwamem2_index = ref_data_types.contains(Constants.RefDataType.BWAMEM2_INDEX) || ref_data_types.contains(Constants.RefDataType.DNA_ALIGNMENT)
        def require_star_index = ref_data_types.contains(Constants.RefDataType.STAR_INDEX) || ref_data_types.contains(Constants.RefDataType.RNA_ALIGNMENT)

        def require_gridss_index = ref_data_types.contains(Constants.RefDataType.GRIDSS_INDEX)
        def require_hmftools_data = ref_data_types.contains(Constants.RefDataType.HMFTOOLS)
        def require_panel_data = ref_data_types.contains(Constants.RefDataType.PANEL)

        if (require_panel_data) {
            if (! params.containsKey('panel') || params.panel == null) {
                require_panel_data = false
                log.warn "Skipping preparing panel specific reference data as --panel CLI argument was not provided"
            } else if (! Constants.PANELS_DEFINED.contains(params.panel.toLowerCase())) {
                require_panel_data = false
                log.warn "Skipping preparing panel specific reference data for custom panel: ${params.panel.toLowerCase()}"
            }
        }

        return [
            prepare_ref_data_only: true,

            require_fasta: require_fasta,
            require_fai: require_fai,
            require_dict: require_dict,
            require_img: require_img,

            require_bwamem2_index: require_bwamem2_index,
            require_star_index: require_star_index,

            require_gridss_index: require_gridss_index,
            require_hmftools_data: require_hmftools_data,
            require_panel_data: require_panel_data,
        ]
    }
}
