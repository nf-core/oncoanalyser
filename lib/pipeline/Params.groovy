package pipeline

import panel.SupportedPanel
import refdata.RefDataType
import refgenome.RefGenomeType
import refgenome.RefGenomeVersion
import refgenome.SupportedGenome
import util.Enums
import util.Messages

class Params {

    public static void parse(Map params) {

        validateRunModes(params)
        validateSequencingType(params)

        validateGenomeAndSetDefaults(params)

        setDefaultHmfData(params)
        validatePanelDataAndSetDefaults(params)
        checkMissingPanelDataPaths(params)

        validateRefDataPathOverrides(params)

        setUmiDefaults(params)
        validateUmiParams(params)
    }

    private static void validateRunModes(Map params) {

        // Pipeline mode
        if (!params.mode) {
            Messages.error(
                "Pipeline mode must be set using the --mode CLI argument or in a configuration file.",
                "",
                "Currently, the available pipeline modes are:",
                Messages.createBulletedList(PipelineMode),
            )
        }

        def pipeline_mode = PipelineMode.fromString((String) params.mode)

        // Purity estimate mode
        if (pipeline_mode == PipelineMode.PURITY_ESTIMATE) {

            if(!params.purity_estimate_mode){
                Messages.error(
                    "A valid purity estimate run mode must be set using the --purity_estimate_mode",
                    "CLI argument or in a configuration file.",
                    "",
                    "Currently, the available modes are:",
                    Messages.createBulletedList(PurityEstimateMode),
                )
            }

            Enums.validateEnumFromString((String) params.purity_estimate_mode, PurityEstimateMode)
        }

        // Prepare reference
        if (pipeline_mode == PipelineMode.PREPARE_REFERENCE && params.ref_data_types == null) {
            Messages.error(
                "CLI argument --ref_data_types is required for mode prepare_reference.",
                "Please specify one or more of the below valid values (separated by commas)",
                Messages.createBulletedList(RefDataType),
            )
        }
    }

    private static void validateGenomeAndSetDefaults(Map params){

        if (!params.genome) {
            Messages.error(
                "Genome must be set using the --genome CLI argument or in a configuration file.",
                "Currently, the supported genomes are: ${params.genomes.keySet().join(", ")}"
            )
        }

        def supported_genome = SupportedGenome.fromString((String) params.genome)
        params.genome_version = (supported_genome != null) ? supported_genome.genomeVersion().getNumericName() : null
        params.genome_type    = (supported_genome != null) ? supported_genome.genomeType() : null

        if(supported_genome)
            return

        if (!params.force_genome) {
            Messages.error(
                "Got unsupported genome: ${params.genome}",
                "",
                "Provide argument --force_genome if you are using a custom genome,",
                "or adjust the --genome argument to one of the supported genomes:",
                Messages.createBulletedList(SupportedGenome)
            )
        }

        if (!params.genome_version) {
            Messages.error(
                "For custom genomes, please provide one of the following values to arg --genome_version:",
                Messages.createBulletedList(RefGenomeVersion.getNumericNames())
            )
        }

        if (!params.genome_type) {
            Messages.error(
                "For custom genomes, please provide of the following values to arg --genome_type:",
                Messages.createBulletedList(RefGenomeType)
            )
        }

        if (params.containsKey('ref_data_genome_alt') && params.ref_data_genome_alt != null) {

            if (params.genome_type != RefGenomeType.ALT) {
                Messages.error("Using a reference genome without ALT contigs but found an .alt file")
            }

            def ref_data_genome_alt_fn = nextflow.Nextflow.file(params.ref_data_genome_alt).getNumericName
            def ref_data_genome_fasta_fn = nextflow.Nextflow.file(params.ref_data_genome_fasta).getNumericName
            if (ref_data_genome_alt_fn != "${ref_data_genome_fasta_fn}.alt") {
                Messages.error(
                    "Found .alt file with filename of ${ref_data_genome_alt_fn} but it is required to match",
                    "reference genome FASTA filename stem: ${ref_data_genome_fasta_fn}.alt"
                )
            }
        }

        params.putIfAbsent('ref_data_genome_alt', null)
        params.putIfAbsent('ref_data_genome_gtf', null)
    }

    private static void validateSequencingType(Map params) {
        Enums.validateEnumFromString((String) params.sequencing_type, SequencingType, false)
    }

    private static void setDefaultHmfData(Map params) {

        if(params.containsKey('ref_data_hmf_data_path'))
            return

        def genome_version = RefGenomeVersion.fromNumericName((String) params.genome_version)
        params.ref_data_hmf_data_path = genome_version.defaultHmfDataPath()
    }

    private static void validatePanelDataAndSetDefaults(Map params){

        def pipeline_mode = PipelineMode.fromString((String) params.mode)
        if (pipeline_mode != PipelineMode.TARGETED)
            return

        if (params.panel == null) {

            Messages.error(
                "A panel is required to be set using the --panel CLI argument or in a ",
                "configuration file when running in targeted mode or panel resource creation mode.",
                "",
                "Currently, panels with built-in support are:",
                Messages.createBulletedList(SupportedPanel.listAll())
            )
        }

        def supported_panel = SupportedPanel.from((String) params.panel, (String) params.genome_version)
        if (supported_panel) {

            if(!params.containsKey('ref_data_panel_data_path')){
                params.ref_data_panel_data_path = supported_panel.defaultDataPath()
            }

        } else {

            if (!params.force_panel) {
                Messages.error(
                    "No built-in support for panel: ${params.panel} ",
                    "",
                    "Provide argument --force_panel if you have a custom panel, or adjust the --panel",
                    "argument to one of the panels configured in the pipeline (case-sensitive):",
                    Messages.createBulletedList(SupportedPanel.listAll())
                )
            }

            if (!params.containsKey('ref_data_panel_data_path')) {
                Messages.error(
                    "If you have a custom panel, provide the directory containing the panel ref data using ",
                    "the CLI argument --ref_data_panel_data_path or in a configuration file"
                )
            }
        }
    }

    private static void checkMissingPanelDataPaths(Map params) {

        def pipeline_mode = PipelineMode.fromString((String) params.mode)
        if (pipeline_mode != PipelineMode.TARGETED)
            return

        def genome_version = RefGenomeVersion.fromNumericName((String) params.genome_version)

        def required_keys = [
            'driver_gene_panel',
            'pon_artefacts',
            'target_region_bed',
            'target_region_normalisation',
        ]

        def optional_keys = [
            'msi_model_error_rates',
            'known_umis',
            'isofox_tpm_norm',
            'isofox_counts',
            'isofox_gc_ratios',
        ]

        for (key in required_keys+optional_keys) {

            def filename = params.panel_data_paths?[params.panel]?[genome_version.getNumericName()]?[key]

            def required_filename_invalid = required_keys.contains(key) && !filename
            def optional_filename_invalid = optional_keys.contains(key) && !filename && filename != []

            if (required_filename_invalid || optional_filename_invalid) {

                def descriptions = []
                descriptions += required_keys.collect { "${it}: Require non-empty path" }
                descriptions += optional_keys.collect { "${it}: Optional non-empty path. If not applicable, set to []" }

                Messages.error(
                    "Panel data filename not defined or misconfigured:",
                    "   params.panel_data_paths.${params.panel}.${genome_version.getNumericName()}.${key} = ${filename}",
                    "",
                    "The below panel data filenames should be configured:",
                    Messages.createBulletedList(descriptions)
                )
            }
        }
    }

    private static void validateRefDataPathOverrides(Map params) {

        def keys = [
            'isofox_counts',
            'isofox_gc_ratios',
            'isofox_gene_ids',
            'isofox_tpm_norm',
            'driver_gene_panel',
            'target_regions_bed',
        ]

        for (key in keys) {
            def path = params[key]

            if (path) {
                nextflow.Nextflow.file(params[key], checkIfExists: true)
            }
        }

    }

    private static void setUmiDefaults(Map params) {

        def umi_type = UmiType.NONE

        if (params.containsKey('umi_type')){

            umi_type = UmiType.fromString((String) params.umi_type)

        } else if (params.containsKey('panel')) {

            def maybe_supported_panel = SupportedPanel.from((String) params.panel, (String) params.genome_version)
            umi_type = UmiType.fromSupportedPanel(maybe_supported_panel)

        }

        // Set defaults if params not set by user
        params.putIfAbsent('fastp_umi_enabled', umi_type.fastpArgs().umiProcessingEnabled())
        params.putIfAbsent('fastp_umi_location', umi_type.fastpArgs().umiLocation())
        params.putIfAbsent('fastp_umi_length', umi_type.fastpArgs().umiLength())
        params.putIfAbsent('fastp_umi_skip', umi_type.fastpArgs().umiSkip())

        params.putIfAbsent('fastq_tools_umi_enabled', umi_type.fastqToolsArgs().umiProcessingEnabled())
        params.putIfAbsent('fastq_tools_umi_delim', umi_type.fastqToolsArgs().umiDelim())

        params.putIfAbsent('redux_umi_enabled', umi_type.reduxArgs().umiProcessingEnabled())
        params.putIfAbsent('redux_umi_duplex_delim', umi_type.reduxArgs().duplexUmiDelim())
    }

    private static void validateUmiParams(Map params) {

        if (params.fastp_umi_enabled && params.fastq_tools_umi_enabled) {
            Messages.error("Either fastp or fastq-tools (not both) can be enabled for UMI processing")
        }

        def fastq_umi_enabled_but_bam_umi_disabled =
            (params.fastp_umi_enabled || params.fastq_tools_umi_enabled) &&
            !params.redux_umi_enabled

        if (fastq_umi_enabled_but_bam_umi_disabled) {
            Messages.error(
                "When FASTQ UMI processing is enabled (with params fastp_umi_enabled or fastp_umi_enabled),",
                "BAM UMI processing with REDUX should also be enabled (with param redux_umi_enabled)"
            )
        }

        // fastp
        def fastp_enabled_but_not_configured =
            params.fastp_umi_enabled &&
            !(params.fastp_umi_location && params.fastp_umi_length && params.fastp_umi_skip >= 0)

        if (fastp_enabled_but_not_configured) {
            Messages.error("fastp UMI processing is enabled but not all of the fastp_umi_* params have been set")
        }

        def fastp_disabled_but_configured =
            !params.fastp_umi_enabled &&
            (params.fastp_umi_location || params.fastp_umi_length || params.fastp_umi_skip >= 0)

        if (fastp_disabled_but_configured) {
            Messages.error("fastp UMI processing is not enabled with param fastp_umi_enabled but detected use of fastp_umi_* params")
        }

        // fastq-tools
        if (params.fastq_tools_umi_enabled && !params.fastq_tools_umi_delim) {
            Messages.error("fastq-tools UMI processing is enabled but param fastq_tools_umi_delim is not set")
        }

        if (!params.fastq_tools_umi_enabled && params.fastq_tools_umi_delim) {
            Messages.error("fastq-tools UMI processing is not enabled with param fastq_tools_umi_enabled but detected use of param fastq_tools_umi_delim")
        }

        // REDUX
        if (!params.redux_umi_enabled && params.redux_umi_duplex_delim) {
            Messages.error("REDUX UMI processing is not enabled with param redux_umi_enabled but detected use of param redux_umi_duplex_delim")
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
