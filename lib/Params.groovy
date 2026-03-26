class Params {

    public static void parseParams(params) {

        validateRunModes(params)

        setGenomeVersionDefaults(params)
        validateGenomeParams(params)

        validatePanelParams(params)

        setUmiDefaults(params)
        validateUmiParams(params)
    }

    private static void error(String... message){

        def message_string = ""
        message_string += "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        message_string += message.join('\n')
        message_string += "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

        throw new RuntimeException(message_string)
    }

    private static String createBulletedList(List<String> items) {
        def bulleted_items = []

        for (item in items) {
            bulleted_items.add(" - ${item}")
        }

        return bulleted_items.join('\n')
    }

    private static void validateRunModes(params) {

        // Pipeline mode
        if (!params.mode) {
            error(
                "Pipeline mode must be set using the --mode CLI argument or in a configuration file.",
                "",
                "Currently, the available pipeline modes are:",
                createBulletedList(Enums.getEnumNames(RunModes.Pipeline)),
            )
        }

        def pipeline_mode = RunModes.Pipeline.fromString(params.mode)

        // Purity estimate mode
        if (pipeline_mode === RunModes.Pipeline.PURITY_ESTIMATE) {

            if(!params.purity_estimate_mode){
                error(
                    "A valid purity estimate run mode must be set using the --purity_estimate_mode",
                    "CLI argument or in a configuration file.",
                    "",
                    "Currently, the available modes are:",
                    createBulletedList(Enums.getEnumNames(RunModes.PurityEstimate)),
                )
            }

            Enums.validateEnumFromString(params.purity_estimate_mode, RunModes.PurityEstimate)
        }

        // Prepare reference
        if (pipeline_mode === RunModes.Pipeline.PREPARE_REFERENCE && params.ref_data_types == null) {
            error(
                "CLI argument --ref_data_types is required for mode prepare_reference.",
                "Please specify one or more of the below valid values (separated by commas)",
                createBulletedList(Enums.getEnumNames(RefData.Type)),
            )
        }

        // Sequencing type
        Enums.validateEnumFromString(params.sequencing_type, RunModes.SequencingType, false)
    }

    private static void setGenomeVersionDefaults(params) {

        if (!params.containsKey('genome_version')) {
            if (RefGenome.GENOMES_VERSION_37.contains(params.genome)) {
                params.genome_version = RefGenome.Version.V37.getNumericName()
            } else if (RefGenome.GENOMES_VERSION_38.contains(params.genome)) {
                params.genome_version = RefGenome.Version.V38.getNumericName()
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
            if (params.genome_version == RefGenome.Version.V37.getNumericName()) {
                params.ref_data_hmf_data_path = RefData.HMF_DATA_37_PATH
            } else if (params.genome_version == RefGenome.Version.V38.getNumericName()) {
                params.ref_data_hmf_data_path = RefData.HMF_DATA_38_PATH
            }
        }

        // Attempt to set default panel data path; make no assumption on valid 'panel' value
        def run_mode = RunModes.Pipeline.fromString(params.mode)
        if ((run_mode === RunModes.Pipeline.TARGETED || run_mode === RunModes.Pipeline.PREPARE_REFERENCE) && params.containsKey('panel')) {
            if (params.panel == 'tso500') {
                if (params.genome_version == RefGenome.Version.V37.getNumericName()) {
                    params.ref_data_panel_data_path = RefData.TSO500_PANEL_37_PATH
                } else if (params.genome_version.toString() == RefGenome.Version.V38.getNumericName()) {
                    params.ref_data_panel_data_path = RefData.TSO500_PANEL_38_PATH
                }
            }
        }
    }

    private static void validateGenomeParams(params) {

        if (!params.ref_data_hmf_data_path) {
            error("CLI argument --ref_data_hmf_data_path must be provided")
        }

        if (!params.genome) {
            error(
                "Genome must be set using the --genome CLI argument or in a configuration file.",
                "Currently, the supported genomes are: ${params.genomes.keySet().join(", ")}"
            )
        }

        if (!RefGenome.GENOMES_SUPPORTED.contains(params.genome) && !params.force_genome) {
            error(
                "Got unsupported genome: ${params.genome}",
                "",
                "Provide argument --force_genome if you are using a custom genome,",
                "or adjust the --genome argument to one of the genomes configured ",
                "in the pipeline: ${RefGenome.GENOMES_SUPPORTED.join(", ")}"
            )
        }

        if (!params.genome_version) {
            error(
                "Genome version wasn't provided and genome '${params.genome}' is not defined in",
                "genome version list. Currently, the list of genomes in the version list include:",
                "${RefGenome.GENOMES_DEFINED.join(", ")}",
            )
        }

        if (!params.genome_type) {
            error(
                "Genome type wasn't provided and genome '${params.genome}' is not defined in",
                "genome type list. Currently, the list of genomes in the type list include:",
                "${RefGenome.GENOMES_DEFINED.join(", ")}",
            )
        }

        if (params.containsKey('ref_data_genome_alt') && params.ref_data_genome_alt !== null) {

            if (params.genome_type != RefGenome.Type.ALT.getName()) {
                error("Using a reference genome without ALT contigs but found an .alt file")
            }

            def ref_data_genome_alt_fn = nextflow.Nextflow.file(params.ref_data_genome_alt).getNumericName
            def ref_data_genome_fasta_fn = nextflow.Nextflow.file(params.ref_data_genome_fasta).getNumericName
            if (ref_data_genome_alt_fn != "${ref_data_genome_fasta_fn}.alt") {
                error(
                    "Found .alt file with filename of ${ref_data_genome_alt_fn} but it is required to match",
                    "reference genome FASTA filename stem: ${ref_data_genome_fasta_fn}.alt"
                )
            }
        }

        if (!params.containsKey('ref_data_genome_alt')) params.ref_data_genome_alt = null
        if (!params.containsKey('ref_data_genome_gtf')) params.ref_data_genome_gtf = null
    }

    private static void validatePanelParams(params) {

        def pipeline_mode = RunModes.Pipeline.fromString(params.mode)

        if (pipeline_mode !== RunModes.Pipeline.TARGETED){
            params.panel = null
            params.ref_data_panel_data_path = null
            return
        }

        def panels_string = createBulletedList(RefData.PANELS_DEFINED)

        if (!params.containsKey('panel') || params.panel === null) {

            error(
                "A panel is required to be set using the --panel CLI argument or in a ",
                "configuration file when running in targeted mode or panel resource creation mode.",
                "",
                "Currently, the supported panels are:",
                panels_string
            )
        }

        if (!RefData.PANELS_DEFINED.contains(params.panel) && (!params.containsKey('force_panel') || !params.force_panel)) {

            error(
                "Got unsupported panel: ${params.panel} ",
                "",
                "Provide argument --force_panel if you have a custom panel,",
                "or adjust the --panel argument to one of the panels configured ",
                "in the pipeline:",
                panels_string
            )
        }
    }

    private static void setUmiDefaults(params) {

        def pipeline_mode = RunModes.Pipeline.fromString(params.mode)

        if (pipeline_mode === RunModes.Pipeline.TARGETED) {

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

    private static void validateUmiParams(params) {

        def fastp_umi_args_set_any = params.fastp_umi_location || params.fastp_umi_length || params.fastp_umi_skip >= 0

        if (fastp_umi_args_set_any && !params.fastp_umi_enabled) {
            error(
                "Detected use of fastp UMI parameters but fastp UMI processing has not been enabled.",
                "Please review your configuration and set the fastp_umi_enabled flag or otherwise ",
                "adjust accordingly."
            )
        }

        def fastp_umi_args_set_all = params.fastp_umi_location && params.fastp_umi_length && params.fastp_umi_skip >= 0

        if (params.fastp_umi_enabled && !fastp_umi_args_set_all) {
            error(
                "Refusing to run fastp UMI processing without having any UMI params configured.",
                "Please review your configuration and appropriately set all fastp_umi_* parameters.",
            )
        }

        if (params.redux_umi_duplex_delim && params.redux_umi_enabled === false) {
            error(
                "Detected use of REDUX UMI parameters but REDUX UMI processing has not been",
                "enabled. Please review your configuration and set the redux_umi_enabled flag or",
                "otherwise adjust accordingly."
            )
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
