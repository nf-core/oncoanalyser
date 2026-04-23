package pipeline

import refdata.RefDataDefaultPaths
import refdata.RefDataType
import refgenome.RefGenomeType
import refgenome.RefGenomeVersion
import refgenome.SupportedGenome
import util.Enums

class Params {

    public static void parse(Map params) {

        validateRunModes(params)
        validateSequencingType(params)

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

    private static void validateRunModes(Map params) {

        // Pipeline mode
        if (!params.mode) {
            error(
                "Pipeline mode must be set using the --mode CLI argument or in a configuration file.",
                "",
                "Currently, the available pipeline modes are:",
                Enums.createBulletedList(PipelineMode),
            )
        }

        def pipeline_mode = PipelineMode.fromString((String) params.mode)

        // Purity estimate mode
        if (pipeline_mode == PipelineMode.PURITY_ESTIMATE) {

            if(!params.purity_estimate_mode){
                error(
                    "A valid purity estimate run mode must be set using the --purity_estimate_mode",
                    "CLI argument or in a configuration file.",
                    "",
                    "Currently, the available modes are:",
                    Enums.createBulletedList(PurityEstimateMode),
                )
            }

            Enums.validateEnumFromString((String) params.purity_estimate_mode, PurityEstimateMode)
        }

        // Prepare reference
        if (pipeline_mode == PipelineMode.PREPARE_REFERENCE && params.ref_data_types == null) {
            error(
                "CLI argument --ref_data_types is required for mode prepare_reference.",
                "Please specify one or more of the below valid values (separated by commas)",
                Enums.createBulletedList(RefDataType),
            )
        }
    }

    private static void validateGenomeAndSetDefaults(Map params){

        if (!params.genome) {
            error(
                "Genome must be set using the --genome CLI argument or in a configuration file.",
                "Currently, the supported genomes are: ${params.genomes.keySet().join(", ")}"
            )
        }

        def supported_genome = SupportedGenome.fromString((String) params.genome)
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
                Enums.createBulletedList(SupportedGenome)
            )
        }

        if (!params.genome_version) {
            error(
                "For custom genomes, please provide one of the following values to arg --genome_version:",
                Enums.createBulletedList(RefGenomeVersion.getNumericNames())
            )
        }

        if (!params.genome_type) {
            error(
                "For custom genomes, please provide of the following values to arg --genome_type:",
                Enums.createBulletedList(RefGenomeType)
            )
        }

        if (params.containsKey('ref_data_genome_alt') && params.ref_data_genome_alt != null) {

            if (params.genome_type != RefGenomeType.ALT) {
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

    private static void validateSequencingType(Map params) {
        Enums.validateEnumFromString((String) params.sequencing_type, SequencingType, false)
    }

    private static void setDefaultHmfData(Map params) {

        if(params.ref_data_hmf_data_path)
            return

        def genome_version = RefGenomeVersion.fromNumericName((String) params.genome_version)
        params.ref_data_hmf_data_path = RefDataDefaultPaths.hmfData(genome_version)
    }

    private static void validatePanelDataAndSetDefaults(Map params){

        def pipeline_mode = PipelineMode.fromString((String) params.mode)

        if (pipeline_mode != PipelineMode.TARGETED && pipeline_mode != PipelineMode.PREPARE_REFERENCE)
            return

        if (params.panel == null) {

            error(
                "A panel is required to be set using the --panel CLI argument or in a ",
                "configuration file when running in targeted mode or panel resource creation mode.",
                "",
                "Currently, panels with built-in support are:",
                Enums.createBulletedList(SupportedPanel)
            )
        }

        def supported_panel = SupportedPanel.fromString((String) params.panel)

        if (!supported_panel) {

            if (params.force_panel)
                return

            error(
                "No built-in supported for panel: ${params.panel} ",
                "",
                "Provide argument --force_panel if you have a custom panel, or adjust the --panel",
                "argument to one of the panels configured in the pipeline (case-sensitive):",
                Enums.createBulletedList(SupportedPanel)
            )
        }

        def genome_version = RefGenomeVersion.fromNumericName((String) params.genome_version)
        if (!supported_panel.hasConfiguredVersion(params, genome_version))
        {
            error("panel(${params.panel}) has no built-in support for refGenomeVersion(${params.genome_version})")
        }

        if (!params.containsKey('ref_data_panel_data_path')) {
            def ref_genome_version = RefGenomeVersion.fromNumericName((String) params.genome_version)
            params.ref_data_panel_data_path = RefDataDefaultPaths.panelData(supported_panel, ref_genome_version)
        }
    }

    private static void setUmiDefaults(Map params) {

        def umi_settings
        if (params.containsKey('panel')) {
            def maybe_supported_panel = SupportedPanel.fromString((String) params.panel)
            umi_settings = UmiSettings.fromSupportedPanel(maybe_supported_panel)
        } else {
            umi_settings = UmiSettings.NONE
        }

        // Set defaults if params not set by user
        params.putIfAbsent('fastp_umi_enabled', umi_settings.fastpArgs().umiProcessingEnabled())
        params.putIfAbsent('fastp_umi_location', umi_settings.fastpArgs().umiLocation())
        params.putIfAbsent('fastp_umi_length', umi_settings.fastpArgs().umiLength())
        params.putIfAbsent('fastp_umi_skip', umi_settings.fastpArgs().umiSkip())

        params.putIfAbsent('fastq_tools_umi_enabled', umi_settings.fastqToolsArgs().umiProcessingEnabled())
        params.putIfAbsent('fastq_tools_umi_delim', umi_settings.fastqToolsArgs().umiDelim())

        params.putIfAbsent('redux_umi_enabled', umi_settings.reduxArgs().umiProcessingEnabled())
        params.putIfAbsent('redux_umi_duplex_delim', umi_settings.reduxArgs().duplexUmiDelim())

        // When fastp UMI is enabled, REDUX UMI should be as well
        params.redux_umi_enabled = params.fastp_umi_enabled
    }

    private static void validateUmiParams(Map params) {

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

        if (params.redux_umi_duplex_delim && params.redux_umi_enabled == false) {
            error(
                "Detected use of REDUX UMI parameters but REDUX UMI processing has not been",
                "enabled. Please review your configuration and set the redux_umi_enabled flag or",
                "otherwise adjust accordingly."
            )
        }
    }

    public static void createStubPlaceholders(Map params) {

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

        params.hmf_data_paths[(String) params.genome_version]
            .each { k, v ->
                fps << "${params.ref_data_hmf_data_path.replaceAll('/$', '')}/${v}"
            }

        if (params.panel != null) {
            params.panel_data_paths[(String) params.panel][(String) params.genome_version]
                .each { k, v ->
                    fps << "${params.ref_data_panel_data_path.replaceAll('/$', '')}/${v}"
                }
        }

        fps.each { fp_str ->
            if (fp_str == null) {
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
