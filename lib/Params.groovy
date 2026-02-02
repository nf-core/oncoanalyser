import nextflow.Nextflow

class Params {

    //
    // Set parameter defaults where required
    //
    public static void setParamsDefaults(params) {

        def default_invalid = false

        // Set defaults common to all run configuration
        if (!params.containsKey('genome_version')) {
            if (Constants.GENOMES_VERSION_37.contains(params.genome)) {
                params.genome_version = '37'
            } else if (Constants.GENOMES_VERSION_38.contains(params.genome)) {
                params.genome_version = '38'
            } else {
                default_invalid = true
            }
        }

        if (!params.containsKey('genome_type')) {
            if (Constants.GENOMES_ALT.contains(params.genome)) {
                params.genome_type = 'alt'
            } else if (Constants.GENOMES_DEFINED.contains(params.genome)) {
                params.genome_type = 'no_alt'
            } else {
                default_invalid = true
            }
        }

        if (!params.containsKey('ref_data_hmf_data_path')) {
            if (params.genome_version.toString() == '37') {
                params.ref_data_hmf_data_path = Constants.HMF_DATA_37_PATH
            } else if (params.genome_version.toString() == '38') {
                params.ref_data_hmf_data_path = Constants.HMF_DATA_38_PATH
            } else {
                default_invalid = true
            }
        }

        // Bad configuration, catch in validateParams
        if (default_invalid) {
            return
        }

        // Set defaults specific to run configuration without attempting to validate

        def run_mode = params.mode ? Enums.getEnumFromString(params.mode, Constants.RunMode) : null
        if (!run_mode) {
            return
        }

        // Attempt to set default panel data path; make no assumption on valid 'panel' value
        if ((run_mode === Constants.RunMode.TARGETED || run_mode === Constants.RunMode.PREPARE_REFERENCE) && params.containsKey('panel')) {
            if (params.panel == 'tso500') {
                if (params.genome_version.toString() == '37') {
                    params.ref_data_panel_data_path = Constants.TSO500_PANEL_37_PATH
                } else if (params.genome_version.toString() == '38') {
                    params.ref_data_panel_data_path = Constants.TSO500_PANEL_38_PATH
                }
            }
        }

        if (run_mode === Constants.RunMode.TARGETED) {

            // When fastp UMI is enabled, REDUX UMI should be as well
            if (params.fastp_umi_enabled && (!params.containsKey('redux_umi_enabled') || !params.redux_umi_enabled)) {
                params.redux_umi_enabled = true
            }

            // Set the REDUX UMI duplex delimiter to '_' when the following conditions are met:
            //   - both fastp and REDUX UMI processing enabled
            //   - fastp is using a duplex UMI location type (per_index or per_read)
            //   - no REDUX duplex delimiter has been set
            def fastp_and_redux_umi_enabled = params.fastp_umi_enabled && params.redux_umi_enabled
            def fastp_duplex_location = params.containsKey('fastp_umi_location') && (params.fastp_umi_location == 'per_index' || params.fastp_umi_location == 'per_read')
            def no_umi_duplex_delim = !params.containsKey('redux_umi_duplex_delim') || !params.redux_umi_duplex_delim
            if (fastp_and_redux_umi_enabled && fastp_duplex_location && no_umi_duplex_delim) {
                params.redux_umi_duplex_delim = '_'
            }

        }

        // Final point to set any default to avoid access to undefined parameters during nf-validation
        if (!params.containsKey('panel')) params.panel = null
        if (!params.containsKey('ref_data_genome_alt')) params.ref_data_genome_alt = null
        if (!params.containsKey('ref_data_genome_gtf')) params.ref_data_genome_gtf = null
        if (!params.containsKey('ref_data_panel_data_path')) params.ref_data_panel_data_path = null

        // Additionally set selected parameters with false-ish truthy values to avoid passing null values as inputs
        if (!params.containsKey('fastp_umi_location')) params.fastp_umi_location = ''
        if (!params.containsKey('fastp_umi_length')) params.fastp_umi_length = 0
        if (!params.containsKey('fastp_umi_skip')) params.fastp_umi_skip = -1
        if (!params.containsKey('redux_umi_duplex_delim')) params.redux_umi_duplex_delim = ''

    }

    public static void validateParams(params, log) {

        if (!params.mode) {
            def run_modes = Enums.getEnumNames(Constants.RunMode).join('\n    - ')
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Run mode must be set using the --mode CLI argument or in a configuration file.\n" +
                "  Currently, the available run modes are:\n" +
                "    - ${run_modes}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        } else {
            Enums.validateEnumFromString(params.mode, Constants.RunMode, log)
        }

        // Genome related

        if (!params.genome) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome must be set using the --genome CLI argument or in a configuration file.\n" +
                "  Currently, the available genome are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        } else if (!params.genomes.containsKey(params.genome)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        if (!Constants.GENOMES_SUPPORTED.contains(params.genome)) {
            if (!params.force_genome) {
                log.error "currently only the GRCh37_hmf and GRCh38_hmf genomes are supported but got ${params.genome}" +
                    ", please adjust the --genome argument accordingly or override with --force_genome."
                Nextflow.exit(1)
            } else {
                log.warn "currently only the GRCh37_hmf and GRCh38_hmf genomes are supported but forcing to " +
                    "proceed with \"${params.genome}\""
            }
        }

        if (!params.genome_version) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome version wasn't provided and genome '${params.genome}' is not defined in   \n" +
                "  genome version list.\n" +
                "  Currently, the list of genomes in the version list include:\n" +
                "  ${Constants.GENOMES_DEFINED.join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        if (!params.genome_type) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome type wasn't provided and genome '${params.genome}' is not defined in      \n" +
                "  genome type list.\n" +
                "  Currently, the list of genomes in the type list include:\n" +
                "  ${Constants.GENOMES_DEFINED.join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        if (!params.ref_data_hmf_data_path) {
            log.error "HMF data path wasn't provided"
            Nextflow.exit(1)
        }

        // Sequencing technology

        Enums.validateEnumFromString(params.sequencing_type, Constants.SequencingType, log, false)

        // UMI parameters

        def fastp_umi_args_set_any = params.fastp_umi_location || params.fastp_umi_length || params.fastp_umi_skip >= 0
        if (fastp_umi_args_set_any && !params.fastp_umi_enabled) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Detected use of fastp UMI parameters but fastp UMI processing has not been enabled.\n" +
                "  Please review your configuration and set the fastp_umi_enabled flag or otherwise " +
                "  adjust accordingly.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        def fastp_umi_args_set_all = params.fastp_umi_location && params.fastp_umi_length && params.fastp_umi_skip >= 0
        if (params.fastp_umi_enabled && !fastp_umi_args_set_all) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Refusing to run fastp UMI processing without having any UMI params configured.\n" +
                "  Please review your configuration and appropriately set all fastp_umi_* parameters.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        if (params.redux_umi_duplex_delim && params.redux_umi_enabled === false) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Detected use of REDUX UMI parameters but REDUX UMI processing has not been\n" +
                "  enabled. Please review your configuration and set the redux_umi_enabled flag or\n" +
                "  otherwise adjust accordingly.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }
    }

    public static getRunConfig(params, inputs, log) {

        def run_mode = Enums.getValidatedEnumFromString(params.mode, Constants.RunMode, log)

        def stages = Processes.getValidatedRunStages(
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

    public static void validateRunSpecificParams(params, run_config, log) {

        // Run mode specific parameters

        if (run_config.run_mode === Constants.RunMode.PREPARE_REFERENCE && params.ref_data_types == null) {

            def ref_data_types = Enums.getEnumNames(Constants.RefDataType).join('\n    - ')
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  CLI argument --ref_data_types is required for mode prepare_reference.\n" +
                "  Please specify one or more of the below valid values (separated by commas)\n" +
                "    - ${ref_data_types}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        if (run_config.run_mode === Constants.RunMode.TARGETED) {

            if (!params.containsKey('panel') || params.panel === null) {

                def panels = Constants.PANELS_DEFINED.join('\n    - ')
                log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  A panel is required to be set using the --panel CLI argument or in a\n" +
                    "  configuration file when running in targeted mode or panel resource creation mode.\n" +
                    "  Currently, the available built-in panels are:\n" +
                    "    - ${panels}\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                Nextflow.exit(1)

            } else if (!Constants.PANELS_DEFINED.contains(params.panel)) {

                if (params.containsKey('force_panel') && params.force_panel) {
                    log.warn "provided panel ${params.panel} does not have built-in support but forcing to proceed"
                } else {
                    def panels = Constants.PANELS_DEFINED.join('\n    - ')
                    log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                        "  The ${params.panel} panel does not have built-in support. Currently, the\n" +
                        "  available supported panels are:\n" +
                        "    - ${panels}\n\n" +
                        "  Please adjust the --panel argument or override with --force_panel.\n" +
                        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                    Nextflow.exit(1)
                }
            }
        }

        if (run_config.run_mode === Constants.RunMode.PURITY_ESTIMATE) {

            if(!params.purity_estimate_mode) {
                def purity_estimate_modes = Enums.getEnumNames(Constants.PurityEstimateRunMode).join('\n    - ')
                log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  A valid purity estimate run mode must be set using the --purity_estimate_mode\n" +
                    "  CLI argument or in a configuration file.\n" +
                    "  Currently, the available run modes are:\n" +
                    "    - ${purity_estimate_modes}\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                Nextflow.exit(1)
            } else {
                Enums.validateEnumFromString(params.purity_estimate_mode, Constants.PurityEstimateRunMode, log)
            }
        }

        if (params.ref_data_genome_alt !== null) {
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

        // NOTE(SW): the following final config checks are performed here since they require additional information
        // regarding processes that are run

        def has_alt_contigs = params.genome_type == 'alt'

        // Ensure that custom genomes with ALT contigs that need indexes built have the required .alt file
        def has_bwa_indexes = (params.ref_data_genome_bwamem2_index && params.ref_data_genome_gridss_index)
        def has_alt_file = params.containsKey('ref_data_genome_alt') && params.ref_data_genome_alt
        def run_bwa_or_gridss_index = run_config.stages.alignment && run_config.has_dna_fastq && !has_bwa_indexes

        if (run_bwa_or_gridss_index && has_alt_contigs && !has_alt_file) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  The genome .alt file is required when building bwa-mem2 or GRIDSS indexes\n" +
                "  for reference genomes containing ALT contigs\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        // Refuse to create STAR index for reference genome containing ALTs, refer to Slack channel
        def run_star_index = run_config.stages.alignment && run_config.has_rna_fastq && !params.ref_data_genome_star_index

        if (run_star_index && has_alt_contigs) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Refusing to create the STAR index for a reference genome with ALT contigs.\n" +
                "  Please review https://github.com/alexdobin/STAR docs or contact us on Slack.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        // Require that an input GTF file is provided when creating STAR index
        if (run_star_index && !params.ref_data_genome_gtf) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Creating a STAR index requires the appropriate genome transcript annotations\n" +
                "  as a GTF file. Please contact us on Slack for further information.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        // Require --isofox_gene_ids argument to be provided in PANEL_RESOURCE_CREATION when RNA inputs are present
        if (run_config.mode === Constants.RunMode.PANEL_RESOURCE_CREATION && run_config.has_rna && !params.isofox_gene_ids) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Running the panel resource creation workflow with RNA requires that the\n" +
                "  --isofox_gene_ids argument is set with an appropriate input file.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }
    }

    public static getPrepConfigFromRunConfig(run_config) {
        return [
            prepare_ref_data_only: false,

            require_fasta: true,
            require_fai: true,
            require_dict: true,
            require_img: true,

            require_bwamem2_index: run_config.has_dna_fastq && run_config.stages.alignment,
            require_star_index: run_config.has_rna_fastq && run_config.stages.alignment,

            require_gridss_index: run_config.has_dna && run_config.mode === Constants.RunMode.WGTS && run_config.stages.virusinterpreter,
            require_hmftools_data: true,
            require_panel_data: run_config.mode === Constants.RunMode.TARGETED,
        ]
    }

    public static getPrepConfigFromParams(params, log) {
        def ref_data_types = params.ref_data_types
            .tokenize(',')
            .collect { Enums.getValidatedEnumFromString(it, Constants.RefDataType, log) }

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
            if (params.panel == null) {
                require_panel_data = false
                log.warn "Skipping preparing panel specific reference data as --panel CLI argument was not provided"
            } else if (!Constants.PANELS_DEFINED.contains(params.panel)) {
                require_panel_data = false
                log.warn "Skipping preparing panel specific reference data for custom panel: ${params.panel}"
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
