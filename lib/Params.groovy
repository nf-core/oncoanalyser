import nextflow.Nextflow

class Params {

    //
    // Set parameter defaults where required
    //
    public static void setDefaults(params) {
        setGenomeVersionDefaults(params)
        setUmiDefaults(params)
        setOtherDefaults(params)
    }

    private static void setGenomeVersionDefaults(params) {

        if (!params.containsKey('genome_version')) {
            if (RefGenome.GENOMES_VERSION_37.contains(params.genome)) {
                params.genome_version = RefGenome.Version.V37.getName()
            } else if (RefGenome.GENOMES_VERSION_38.contains(params.genome)) {
                params.genome_version = RefGenome.Version.V38.getName()
            }
        }

        if (!params.containsKey('genome_type')) {
            if (RefGenome.GENOMES_ALT.contains(params.genome)) {
                params.genome_type = RefGenome.Type.ALT.getName()
            } else if (RefGenome.GENOMES_DEFINED.contains(params.genome)) {
                params.genome_type = RefGenome.Type.NO_ALT.getName()
            }
        }

        if (!params.containsKey('ref_data_hmf_data_path')) {
            if (params.genome_version == RefGenome.Version.V37.getName()) {
                params.ref_data_hmf_data_path = RefData.HMF_DATA_37_PATH
            } else if (params.genome_version == RefGenome.Version.V38.getName()) {
                params.ref_data_hmf_data_path = RefData.HMF_DATA_38_PATH
            }
        }

        // Attempt to set default panel data path; make no assumption on valid 'panel' value
        def run_mode = RunModes.Pipeline.fromString(params.mode)
        if ((run_mode === RunModes.Pipeline.TARGETED || run_mode === RunModes.Pipeline.PREPARE_REFERENCE) && params.containsKey('panel')) {
            if (params.panel == 'tso500') {
                if (params.genome_version == RefGenome.Version.V37.getName()) {
                    params.ref_data_panel_data_path = RefData.TSO500_PANEL_37_PATH
                } else if (params.genome_version.toString() == RefGenome.Version.V38.getName()) {
                    params.ref_data_panel_data_path = RefData.TSO500_PANEL_38_PATH
                }
            }
        }
    }

    private static void setUmiDefaults(params) {

        def run_mode = RunModes.Pipeline.fromString(params.mode)

        if (run_mode === RunModes.Pipeline.TARGETED) {

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

        // Additionally set selected parameters with false-ish truthy values to avoid passing null values as inputs
        if (!params.containsKey('fastp_umi_location')) params.fastp_umi_location = ''
        if (!params.containsKey('fastp_umi_length')) params.fastp_umi_length = 0
        if (!params.containsKey('fastp_umi_skip')) params.fastp_umi_skip = -1
        if (!params.containsKey('redux_umi_duplex_delim')) params.redux_umi_duplex_delim = ''
    }

    private static void setOtherDefaults(params) {
        // Final point to set any default to avoid access to undefined parameters during nf-validation
        if (!params.containsKey('panel')) params.panel = null
        if (!params.containsKey('ref_data_genome_alt')) params.ref_data_genome_alt = null
        if (!params.containsKey('ref_data_genome_gtf')) params.ref_data_genome_gtf = null
        if (!params.containsKey('ref_data_panel_data_path')) params.ref_data_panel_data_path = null
    }

    public static void validateParams(params, log) {

        if (!params.mode) {
            def run_modes = Enums.getEnumNames(RunModes.Pipeline).join('\n    - ')
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Run mode must be set using the --mode CLI argument or in a configuration file.\n" +
                "  Currently, the available run modes are:\n" +
                "    - ${run_modes}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        } else {
            Enums.validateEnumFromString(params.mode, RunModes.Pipeline)
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

        if (!RefGenome.GENOMES_SUPPORTED.contains(params.genome)) {
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
                "  ${RefGenome.GENOMES_DEFINED.join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        if (!params.genome_type) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome type wasn't provided and genome '${params.genome}' is not defined in      \n" +
                "  genome type list.\n" +
                "  Currently, the list of genomes in the type list include:\n" +
                "  ${RefGenome.GENOMES_DEFINED.join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        if (!params.ref_data_hmf_data_path) {
            log.error "HMF data path wasn't provided"
            Nextflow.exit(1)
        }

        // Sequencing technology

        Enums.validateEnumFromString(params.sequencing_type, RunModes.SequencingType, false)

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

    public static getRunConfig(params, log) {

        def run_mode = RunModes.Pipeline.fromString(params.mode)

        def stages = RunStage.getValidatedRunStages(
            params.processes_include,
            params.processes_exclude,
            params.processes_manual,
            log,
        )

        return [
            mode: run_mode,
            stages: stages,
        ]
    }

    public static void validateRunSpecificParams(params, run_config, log) {

        // Run mode specific parameters

        if (run_config.run_mode === RunModes.Pipeline.PREPARE_REFERENCE && params.ref_data_types == null) {

            def ref_data_types = Enums.getEnumNames(RefData.Type).join('\n    - ')
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  CLI argument --ref_data_types is required for mode prepare_reference.\n" +
                "  Please specify one or more of the below valid values (separated by commas)\n" +
                "    - ${ref_data_types}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }

        if (run_config.run_mode === RunModes.Pipeline.TARGETED) {

            if (!params.containsKey('panel') || params.panel === null) {

                def panels = RefData.PANELS_DEFINED.join('\n    - ')
                log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  A panel is required to be set using the --panel CLI argument or in a\n" +
                    "  configuration file when running in targeted mode or panel resource creation mode.\n" +
                    "  Currently, the available built-in panels are:\n" +
                    "    - ${panels}\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                Nextflow.exit(1)

            } else if (!RefData.PANELS_DEFINED.contains(params.panel)) {

                if (params.containsKey('force_panel') && params.force_panel) {
                    log.warn "provided panel ${params.panel} does not have built-in support but forcing to proceed"
                } else {
                    def panels = RefData.PANELS_DEFINED.join('\n    - ')
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

        if (run_config.run_mode === RunModes.Pipeline.PURITY_ESTIMATE) {

            if(!params.purity_estimate_mode) {
                def purity_estimate_modes = Enums.getEnumNames(RunModes.PurityEstimate).join('\n    - ')
                log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  A valid purity estimate run mode must be set using the --purity_estimate_mode\n" +
                    "  CLI argument or in a configuration file.\n" +
                    "  Currently, the available run modes are:\n" +
                    "    - ${purity_estimate_modes}\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                Nextflow.exit(1)
            } else {
                Enums.validateEnumFromString(params.purity_estimate_mode, RunModes.PurityEstimate)
            }
        }

        if (params.ref_data_genome_alt !== null) {
            if (params.genome_type != RefGenome.Type.ALT.getName()) {
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

        def has_alt_contigs = params.genome_type == RefGenome.Type.ALT.getName()

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
        if (run_config.mode === RunModes.Pipeline.PANEL_RESOURCE_CREATION && run_config.has_rna && !params.isofox_gene_ids) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Running the panel resource creation workflow with RNA requires that the\n" +
                "  --isofox_gene_ids argument is set with an appropriate input file.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.exit(1)
        }
    }

    public static void createStubPlaceholders(params) {

        def fps = [
            params.ref_data_genome_alt,
            params.ref_data_genome_bwamem2_index,
            params.ref_data_genome_dict,
            params.ref_data_genome_fai,
            params.ref_data_genome_fasta,
            params.ref_data_genome_gridss_index,
            params.ref_data_genome_gtf,
            params.ref_data_genome_star_index,
        ]

        params.hmf_data_paths[params.genome_version.toString()]
            .each { k, v ->
                fps << "${params.ref_data_hmf_data_path.replaceAll('/$', '')}/${v}"
            }

        if (params.panel !== null) {
            params.panel_data_paths[params.panel][params.genome_version.toString()]
                .each { k, v ->
                    fps << "${params.ref_data_panel_data_path.replaceAll('/$', '')}/${v}"
                }
        }

        fps.each { fp_str ->
            if (fp_str === null) {
                return
            }

            def fp = fp_str ? nextflow.Nextflow.file(fp_str) : []

            if (!fp_str || fp.exists()) {
                return
            }

            if (fp_str.endsWith('/')) {
                fp.mkdirs()
            } else {
                fp.getParent().mkdirs()
                fp.toFile().createNewFile()
            }
        }
    }
}
