class Params {

    public static void parseParams(params) {

        validateRunModes(params)

        validateGenomeAndSetDefaults(params)

        setDefaultHmfData(params)
        validatePanelDataAndSetDefaults(params)

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

    private static void validateGenomeAndSetDefaults(params) {
        if (!params.genome) {
            error(
                "Genome must be set using the --genome CLI argument or in a configuration file.",
                "Currently, the supported genomes are: ${params.genomes.keySet().join(", ")}"
            )
        }

        def supported_genome = RefGenome.SupportedGenome.fromString(params.genome)
        params.genome_version = (supported_genome != null) ? supported_genome.getVersion().getNumericName() : null
        params.genome_type    = (supported_genome != null) ? supported_genome.getType() : null

        if(supported_genome)
            return

        if (!params.force_genome) {
            error(
                "Got unsupported genome: ${params.genome}",
                "",
                "Provide argument --force_genome if you are using a custom genome,",
                "or adjust the --genome argument to one of the supported genomes:",
                createBulletedList(RefGenome.SupportedGenome.getNames())
            )
        }

        if (!params.genome_version) {
            error(
                "For custom genomes, please provide one of the following values to arg --genome_version:",
                createBulletedList(RefGenome.Version.getNumericNames())
            )
        }

        if (!params.genome_type) {
            error(
                "For custom genomes, please provide of the following values to arg --genome_type:",
                createBulletedList(RefGenome.Type.getNames())
            )
        }

        if (params.containsKey('ref_data_genome_alt') && params.ref_data_genome_alt != null) {

            if (params.genome_type != RefGenome.Type.ALT) {
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

    private static void setDefaultHmfData(params) {

        if(params.ref_data_hmf_data_path)
            return

        def genome_version = RefGenome.Version.fromNumericName(params.genome_version)
        params.ref_data_hmf_data_path = RefData.getDefaultHmfDataPath(genome_version)
    }

    private static void validatePanelDataAndSetDefaults(params){

        def pipeline_mode = RunModes.Pipeline.fromString(params.mode)

        if (pipeline_mode != RunModes.Pipeline.TARGETED && pipeline_mode != RunModes.Pipeline.PREPARE_REFERENCE)
            return

        if (params.panel == null) {

            error(
                "A panel is required to be set using the --panel CLI argument or in a ",
                "configuration file when running in targeted mode or panel resource creation mode.",
                "",
                "Currently, the supported panels are:",
                createBulletedList(RefData.SupportedPanel.getNames())
            )
        }

        def supported_panel = RefData.SupportedPanel.fromString(params.panel)

        if (supported_panel && !params.containsKey('ref_data_panel_data_path')) {
            def ref_genome_version = RefGenome.Version.fromNumericName(params.genome_version)
            params.ref_data_panel_data_path = RefData.getDefaultPanelDataPath(supported_panel, ref_genome_version)
        }

        if (!supported_panel && !params.force_panel) {

            error(
                "Got unsupported panel: ${params.panel} ",
                "",
                "Provide argument --force_panel if you have a custom panel,",
                "or adjust the --panel argument to one of the panels configured ",
                "in the pipeline:",
                createBulletedList(RefData.SupportedPanel.getNames())
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
