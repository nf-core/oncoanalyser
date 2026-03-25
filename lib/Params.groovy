import nextflow.Nextflow

class Params {

    public static void parseParams(params) {

        validateRunModes(params)

        setGenomeVersionDefaults(params)
        validateGenomeParams(params)

        validatePanelParams(params)

        setUmiDefaults(params)
        validateUmiParams(params)
    }

    private static void validateRunModes(params) {

        // Pipeline mode
        if (!params.mode) {

            def pipeline_modes = Enums.getEnumNames(RunModes.Pipeline).join('\n    - ')
            throw new IllegalStateException(
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Pipeline mode must be set using the --mode CLI argument or in a configuration file.\n" +
                "  Currently, the available pipeline modes are:\n" +
                "    - ${pipeline_modes}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            )
        }

        def pipeline_mode = RunModes.Pipeline.fromString(params.mode)

        // Purity estimate mode
        if (pipeline_mode === RunModes.Pipeline.PURITY_ESTIMATE) {

            if(!params.purity_estimate_mode){
                def purity_estimate_modes = Enums.getEnumNames(RunModes.PurityEstimate).join('\n    - ')
                throw new IllegalStateException(
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  A valid purity estimate run mode must be set using the --purity_estimate_mode\n" +
                    "  CLI argument or in a configuration file.\n" +
                    "  Currently, the available run modes are:\n" +
                    "    - ${purity_estimate_modes}\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                )
            }

            Enums.validateEnumFromString(params.purity_estimate_mode, RunModes.PurityEstimate)
        }

        // Prepare reference
        if (pipeline_mode === RunModes.Pipeline.PREPARE_REFERENCE && params.ref_data_types == null) {

            def ref_data_types = Enums.getEnumNames(RefData.Type).join('\n    - ')
            throw new IllegalStateException(
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  CLI argument --ref_data_types is required for mode prepare_reference.\n" +
                "  Please specify one or more of the below valid values (separated by commas)\n" +
                "    - ${ref_data_types}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            )
        }

        // Sequencing type
        Enums.validateEnumFromString(params.sequencing_type, RunModes.SequencingType, false)

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

    private static void validateGenomeParams(params) {

        if (!params.ref_data_hmf_data_path) {
            throw new IllegalArgumentException("CLI argument --ref_data_hmf_data_path must be provided")
        }

        //
        if (!params.genome) {
            throw new IllegalStateException(
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome must be set using the --genome CLI argument or in a configuration file.\n" +
                "  Currently, the available genome are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            )

        } else if (!params.genomes.containsKey(params.genome)) {
            throw new IllegalStateException(
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            )
        }

        if (!RefGenome.GENOMES_SUPPORTED.contains(params.genome) && !params.force_genome) {
            throw new IllegalStateException(
                "currently only the GRCh37_hmf and GRCh38_hmf genomes are supported but got ${params.genome}" +
                ", please adjust the --genome argument accordingly or override with --force_genome."
            )
        }

        if (!params.genome_version) {
            throw new IllegalStateException(
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome version wasn't provided and genome '${params.genome}' is not defined in   \n" +
                "  genome version list.\n" +
                "  Currently, the list of genomes in the version list include:\n" +
                "  ${RefGenome.GENOMES_DEFINED.join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            )
        }

        // Alt genomes
        if (!params.genome_type) {
            throw new IllegalStateException(
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome type wasn't provided and genome '${params.genome}' is not defined in\n" +
                "  genome type list.\n" +
                "  Currently, the list of genomes in the type list include:\n" +
                "  ${RefGenome.GENOMES_DEFINED.join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            )
        }

        if (params.containsKey('ref_data_genome_alt') && params.ref_data_genome_alt !== null) {

            if (params.genome_type != RefGenome.Type.ALT.getName()) {
                throw new IllegalStateException("Using a reference genome without ALT contigs but found an .alt file")
            }

            def ref_data_genome_alt_fn = nextflow.Nextflow.file(params.ref_data_genome_alt).name
            def ref_data_genome_fasta_fn = nextflow.Nextflow.file(params.ref_data_genome_fasta).name
            if (ref_data_genome_alt_fn != "${ref_data_genome_fasta_fn}.alt") {
                throw new IllegalStateException(
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  Found .alt file with filename of ${ref_data_genome_alt_fn} but it is required to match\n" +
                    "  reference genome FASTA filename stem: ${ref_data_genome_fasta_fn}.alt\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
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

        def panels_string = RefData.PANELS_DEFINED.join('\n    - ')

        if (!params.containsKey('panel') || params.panel === null) {

            throw new IllegalStateException(
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  A panel is required to be set using the --panel CLI argument or in a\n" +
                "  configuration file when running in targeted mode or panel resource creation mode.\n" +
                "  Currently, the available built-in panels are:\n" +
                "    - ${panels_string}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            )
        }

        if (!RefData.PANELS_DEFINED.contains(params.panel) && (!params.containsKey('force_panel') || !params.force_panel)) {

            throw new IllegalStateException(
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  The ${params.panel} panel does not have built-in support. Currently, the\n" +
                "  available supported panels are:\n" +
                "    - ${panels_string}\n\n" +
                "  Please adjust the --panel argument or override with --force_panel.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
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
            throw new IllegalStateException(
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Detected use of fastp UMI parameters but fastp UMI processing has not been enabled.\n" +
                "  Please review your configuration and set the fastp_umi_enabled flag or otherwise " +
                "  adjust accordingly.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            )
        }

        def fastp_umi_args_set_all = params.fastp_umi_location && params.fastp_umi_length && params.fastp_umi_skip >= 0

        if (params.fastp_umi_enabled && !fastp_umi_args_set_all) {
            throw new IllegalStateException(
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Refusing to run fastp UMI processing without having any UMI params configured.\n" +
                "  Please review your configuration and appropriately set all fastp_umi_* parameters.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            )
        }

        if (params.redux_umi_duplex_delim && params.redux_umi_enabled === false) {
            throw new IllegalStateException(
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Detected use of REDUX UMI parameters but REDUX UMI processing has not been\n" +
                "  enabled. Please review your configuration and set the redux_umi_enabled flag or\n" +
                "  otherwise adjust accordingly.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
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
