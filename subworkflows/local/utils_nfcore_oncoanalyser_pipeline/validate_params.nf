//
// Parameter default and validation helpers for the nf-core/oncoanalyser pipeline
//

include { getSamples               } from './accessors'
include { getTumorDnaSample        } from './accessors'
include { hasAlignmentInput        } from './accessors'
include { hasInput                 } from './accessors'
include { hasNormalDna             } from './accessors'
include { hasNormalDnaAlignment    } from './accessors'
include { hasNormalDnaFastq        } from './accessors'
include { hasTumorDna              } from './accessors'
include { hasTumorDnaFastq         } from './accessors'
include { hasTumorRna              } from './accessors'
include { hasTumorRnaFastq         } from './accessors'
include { getRunStages             } from './processes'
include { FileType                 } from './types'
include { RefDataType              } from './types'
include { RunMode                  } from './types'
include { UmiType                  } from './types'
include { getEnumFromString        } from './utils'
include { getEnumFromStringOrFail  } from './utils'
include { getEnumNames             } from './utils'
include { getRunMode               } from './utils'

//
// Set parameter defaults where required
//
def setParamsDefaults(params, log) {

    def default_invalid = setCommonDefaults(params)

    // Bad configuration, catch in validateParams
    if (default_invalid) {
        return
    }

    // Set defaults specific to run configuration without attempting to validate
    def run_mode
    if (params.mode != null) {
        run_mode = getRunMode(params.mode, log)
    } else {
        // Bad configuration, catch in validateParams
        return
    }

    setRunModeDefaults(params, run_mode)

    def stages = getRunStages(
        params.processes_include,
        params.processes_exclude,
        params.processes_manual,
        log,
    )

    setUmiDefaults(params, log)

    setParamPlaceholderDefaults(params)
}

def setCommonDefaults(params) {

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

    return default_invalid
}

def setRunModeDefaults(params, run_mode) {

    // Attempt to set default panel data path; make no assumption on valid 'panel' value
    if (run_mode == RunMode.TARGETED || run_mode == RunMode.PREPARE_REFERENCE) {

        if (params.panel != null) {

            if (params.panel.toLowerCase() == 'tso500') {
                if (params.genome_version.toString() == '37') {
                    params.ref_data_panel_data_path = "${params.ref_data_base}/${Constants.TSO500_PANEL_37_PATH}"
                }
            }

        }

    }
}

def setUmiDefaults(params, log) {

    //
    // Resolve UMI type and set UMI parameters
    //
    def umi_type
    if (params.containsKey('umi_type') && params.umi_type) {
        umi_type = getEnumFromString(params.umi_type, UmiType)
    } else if (params.panel != null && Constants.PANELS_DEFINED.contains(params.panel.toLowerCase())) {
        if (params.panel.toLowerCase() == 'tso500') {
            umi_type = UmiType.TSO500
        }
    }

    if (umi_type == UmiType.KAPA) {
          params.fastp_umi_enabled = true
          params.fastp_umi_location = 'per_read'
          params.fastp_umi_length = 6
          params.fastp_umi_skip = 0
          params.redux_umi_enabled = true
          params.redux_umi_duplex_delim = '_'
    } else if (umi_type == UmiType.MSK) {
          params.fastq_tools_umi_enabled = true
          params.fastq_tools_umi_delim = '+'
          params.redux_umi_enabled = true
          params.redux_umi_duplex_delim = ''
    } else if (umi_type == UmiType.TSO500) {
          params.redux_umi_enabled = true
          params.redux_umi_duplex_delim = '+'
    } else if (umi_type == UmiType.TWIST) {
          params.fastp_umi_enabled = true
          params.fastp_umi_location = 'per_read'
          params.fastp_umi_length = 7
          params.fastp_umi_skip = 0
          params.redux_umi_enabled = true
          params.redux_umi_duplex_delim = '_'
    }
}

def setParamPlaceholderDefaults(params) {

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
def validateParams(params, log) {

    validateGenomeParams(params, log)

    def run_mode = validateRunModeParams(params, log)

    validatePrepareReferenceParams(params, run_mode, log)

    validateWgtsParams(params, run_mode, log)

    validateTargetedParams(params, run_mode, log)

    validatePurityEstimateParams(params, run_mode, log)

    validateAltContigsParams(params, log)

    validateUmiParams(params, log)
}

def validateGenomeParams(params, log) {

    // Common parameters

    if (! params.genome) {
        log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome must be set using the --genome CLI argument or in a configuration file.\n" +
            "  Currently, the available genomes are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    } else if (! params.genomes.containsKey(params.genome)) {
        log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genomes are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

    if (! Constants.GENOMES_SUPPORTED.contains(params.genome)) {
        if (! params.force_genome) {
            log.error "currently only the GRCh37_hmf and GRCh38_hmf genomes are supported but got ${params.genome}" +
                ", please adjust the --genome argument accordingly or override with --force_genome."
            exit 1
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
        exit 1
    }

    if (! params.genome_type) {
        log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome type wasn't provided and genome '${params.genome}' is not defined in      \n" +
            "  genome type list.\n" +
            "  Currently, the list of genomes in the type list include:\n" +
            "  ${Constants.GENOMES_DEFINED.join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

    if (! params.ref_data_hmf_data_path) {
        log.error "HMF data path wasn't provided"
        exit 1
    }

    if (! params.hmf_data_paths.containsKey(params.genome_version.toString())) {
        def hmf_data_versions = params.hmf_data_paths.keySet().join('\n    - ')
        log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Did not find path definitions for the provided 'hmf_data_paths' and\n" +
            "  genome version ${params.genome_version}. Please check your configuration.\n" +
            "  Found the following genome version definitions for 'hmf_data_paths':\n" +
            "    - ${hmf_data_versions}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }
}

def validateRunModeParams(params, log) {

    // Run configuration specific parameters

    if (! params.mode) {
        def run_modes = getEnumNames(RunMode).join('\n    - ')
        log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Run mode must be set using the --mode CLI argument or in a configuration  \n" +
            "  file.\n" +
            "  Currently, the available run modes are:\n" +
            "    - ${run_modes}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

    return getRunMode(params.mode, log)
}

def validatePrepareReferenceParams(params, run_mode, log) {

    if (run_mode == RunMode.PREPARE_REFERENCE && params.ref_data_types == null) {

        def ref_data_types = getEnumNames(RefDataType).join('\n    - ')

        log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  CLI argument --ref_data_types is required for mode prepare_reference.\n" +
            "  Please specify one or more of the below valid values (separated by commas)\n" +
            "    - ${ref_data_types}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }
}

def validateWgtsParams(params, run_mode, log) {

    if (run_mode == RunMode.WGTS) {

        // We allow fastq-tools UMI processing in WGTS but this requires that the user manually set the 'known_umis' file in hmf_data, enforce here
        def has_known_umis = params.hmf_data_paths[params.genome_version.toString()].containsKey('known_umis')
        if (params.fastq_tools_umi_enabled && ! has_known_umis) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  The fastq-tools process is enabled but the required 'known_umis' data path was\n" +
                "  not configured in the respective 'hmf_data_paths' for genome version ${params.genome_version.toString()}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            exit 1
        }

    }
}

def validateTargetedParams(params, run_mode, log) {

    if (run_mode == RunMode.TARGETED) {

        if (! params.containsKey('panel') || params.panel == null) {

            def panels = Constants.PANELS_DEFINED.join('\n    - ')
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  A panel is required to be set using the --panel CLI argument or in a\n" +
                "  configuration file when running in targeted mode or panel resource creation mode.\n" +
                "  Currently, the available built-in panels are:\n" +
                "    - ${panels}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            exit 1

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
                exit 1
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
            exit 1
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
            exit 1
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
            exit 1
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
            exit 1
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
            exit 1
        }

        // Require known_umis in panel data when fastq-tools UMI processing is enabled
        if (params.fastq_tools_umi_enabled && ! panel_data_paths['known_umis']) {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  The fastq-tools process is enabled but the required 'known_umis' data path was\n" +
                "  not configured in the panel data paths for panel ${params.panel.toLowerCase()}\n" +
                "  (${params.genome_version}).\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            exit 1
        }

    }
}

def validatePurityEstimateParams(params, run_mode, log) {

    if (run_mode == RunMode.PURITY_ESTIMATE) {

        def purity_estimate_modes = [RunMode.WGTS, RunMode.TARGETED]

        def purity_mode_enum = ! params.purity_estimate_mode ? null : getEnumFromString(params.purity_estimate_mode, RunMode)

        if (! purity_mode_enum || ! purity_estimate_modes.contains(purity_mode_enum)) {

            def purity_estimate_modes_str = purity_estimate_modes.collect { e -> e.name().toLowerCase() }.join('\n    - ')

            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  A valid purity estimate run mode must be set using the --purity_estimate_mode\n" +
                "  CLI argument or in a configuration file.\n" +
                "  Currently, the available run modes are:\n" +
                "    - ${purity_estimate_modes_str}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            exit 1
        }
    }
}

def validateAltContigsParams(params, log) {

    if (params.ref_data_genome_alt != null) {
        if (params.genome_type != 'alt') {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Using a reference genome without ALT contigs but found an .alt file\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            exit 1
        }

        def ref_data_genome_alt_fn = nextflow.Nextflow.file(params.ref_data_genome_alt).name
        def ref_data_genome_fasta_fn = nextflow.Nextflow.file(params.ref_data_genome_fasta).name
        if (ref_data_genome_alt_fn != "${ref_data_genome_fasta_fn}.alt") {
            log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Found .alt file with filename of ${ref_data_genome_alt_fn} but it is required to match\n" +
                "  reference genome FASTA filename stem: ${ref_data_genome_fasta_fn}.alt\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            exit 1
        }

    }
}

def validateUmiParams(params, log) {

    // UMI parameters
    if (params.fastp_umi_enabled && params.fastq_tools_umi_enabled) {
        log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
              "  UMI processing with either fastp or fastq-tools but not both can be enabled\n" +
              "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

    if ((params.fastp_umi_enabled || params.fastq_tools_umi_enabled) && ! params.redux_umi_enabled) {
        log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
              "  When FASTQ UMI processing is enabled (via fastp_umi_enabled or fastq_tools_umi_enabled),\n" +
              "  REDUX UMI processing must also be enabled with redux_umi_enabled\n" +
              "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

    def fastp_umi_args_set_any = params.fastp_umi_location || params.fastp_umi_length || params.fastp_umi_skip >= 0
    if (fastp_umi_args_set_any && ! params.fastp_umi_enabled) {
        log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Detected use of fastp UMI parameters but fastp UMI processing has not been enabled.\n" +
            "  Please review your configuration and set the fastp_umi_enabled flag or otherwise " +
            "  adjust accordingly.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

    def fastp_umi_args_set_all = params.fastp_umi_location && params.fastp_umi_length && params.fastp_umi_skip >= 0
    if (params.fastp_umi_enabled && ! fastp_umi_args_set_all) {
        log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Refusing to run fastp UMI processing without having all UMI params configured.\n" +
            "  Please review your configuration and appropriately set all fastp_umi_* parameters.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

    if (params.fastq_tools_umi_delim && ! params.fastq_tools_umi_enabled) {
        log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Detected use of fastq-tools UMI parameter fastq_tools_umi_delim but fastq-tools UMI\n" +
            "  processing has not been enabled. Please review your configuration and set the\n" +
            "  fastq_tools_umi_enabled or otherwise adjust accordingly.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

    if (params.fastq_tools_umi_enabled && ! params.fastq_tools_umi_delim) {
        log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Refusing to run fastq-tools UMI processing without fastq_tools_umi_delim configured.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

    if (params.redux_umi_duplex_delim && params.redux_umi_enabled == false) {
        log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Detected use of REDUX UMI parameters but REDUX UMI processing has not been\n" +
            "  enabled. Please review your configuration and set the redux_umi_enabled flag or\n" +
            "  otherwise adjust accordingly.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

}

def getRunConfig(params, inputs, log) {

    def run_mode = getRunMode(params.mode, log)

    def stages = getRunStages(
        params.processes_include,
        params.processes_exclude,
        params.processes_manual,
        log,
    )

    return [
        mode: run_mode,
        stages: stages,
        has_dna: inputs.any { hasTumorDna(it) },
        has_rna: inputs.any { hasTumorRna(it) },
        has_rna_fastq: inputs.any { hasTumorRnaFastq(it) },
        has_dna_fastq: inputs.any { hasTumorDnaFastq(it) || hasNormalDnaFastq(it) },
    ]
}

def getPrepConfigFromSamplesheet(run_config) {
    return [
        prepare_ref_data_only: false,

        require_fasta: true,
        require_fai: true,
        require_dict: true,
        require_img: true,

        require_bwamem2_index: run_config.has_dna_fastq && run_config.stages.alignment,
        require_star_index: run_config.has_rna_fastq && run_config.stages.alignment,

        require_gridss_index: run_config.has_dna && run_config.mode == RunMode.WGTS && run_config.stages.virusinterpreter,
        require_hmftools_data: true,
        require_panel_data: run_config.mode == RunMode.TARGETED,
    ]
}

def getPrepConfigFromCli(params, log) {
    def ref_data_types = params.ref_data_types.tokenize(',').collect {
            return getEnumFromStringOrFail(it, RefDataType, 'ref data type', log)
        }

    if (
        ref_data_types.contains(RefDataType.WGS) ||
        ref_data_types.contains(RefDataType.WTS) ||
        ref_data_types.contains(RefDataType.TARGETED)
    ) {
        ref_data_types += [
            RefDataType.FASTA,
            RefDataType.FAI,
            RefDataType.DICT,
            RefDataType.IMG,
            RefDataType.HMFTOOLS
        ]
    }

    if (ref_data_types.contains(RefDataType.WGS)) {
        ref_data_types += [RefDataType.GRIDSS_INDEX]
    }

    if (ref_data_types.contains(RefDataType.TARGETED)) {
        ref_data_types += [RefDataType.PANEL]
    }

    def require_fasta = ref_data_types.contains(RefDataType.FASTA)
    def require_fai = ref_data_types.contains(RefDataType.FAI)
    def require_dict = ref_data_types.contains(RefDataType.DICT)
    def require_img = ref_data_types.contains(RefDataType.IMG)

    def require_bwamem2_index = ref_data_types.contains(RefDataType.BWAMEM2_INDEX) || ref_data_types.contains(RefDataType.DNA_ALIGNMENT)
    def require_star_index = ref_data_types.contains(RefDataType.STAR_INDEX) || ref_data_types.contains(RefDataType.RNA_ALIGNMENT)

    def require_gridss_index = ref_data_types.contains(RefDataType.GRIDSS_INDEX)
    def require_hmftools_data = ref_data_types.contains(RefDataType.HMFTOOLS)
    def require_panel_data = ref_data_types.contains(RefDataType.PANEL)

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

def validateInput(inputs, run_config, params, log) {

    inputs.each { case_record ->

        // Require ALN or ALN_REDUX or REDUX_DIR or FASTQs for each defined sample
        getSamples(case_record).each { sample ->
            if (! hasAlignmentInput(sample)) {
                log.error "no alignments (or REDUX alignments / directory) nor FASTQ files provided for ${case_record.case_id} ${sample.sample_type}/${sample.sequence_type}\n\n" +
                    "NB: At least one of these files is required as they are the basis to determine input sample type."
                exit 1
            }
        }

        // Do not allow donor sample without normal sample
        if (case_record.donor_dna_samples && ! case_record.normal_dna_samples) {
            log.error "a donor sample but no normal sample was found for ${case_record.case_id}\n\n" +
                "Analysis with a donor sample requires a normal sample."
            exit 1
        }

        // Longitudinal samples require a primary normal DNA alignment when AMBER input is provided
        if (case_record.longitudinal_samples && hasInput(getTumorDnaSample(case_record), FileType.AMBER_DIR) && ! hasNormalDnaAlignment(case_record)) {
            log.error "AMBER input was provided without the required primary normal DNA BAM for ${case_record.case_id}"
            exit 1
        }

        // Apply some required restrictions to targeted mode
        if (run_config.mode == RunMode.TARGETED) {

            // Do not allow donor DNA
            if (case_record.donor_dna_samples) {
                log.error "targeted mode is not compatible with the donor DNA BAM/CRAM provided for ${case_record.case_id}\n\n" +
                    "The targeted workflow supports only tumor and normal DNA BAM/CRAMs (and tumor RNA BAM/CRAMs for TSO500)"
                exit 1
            }

            // Do not allow only tumor RNA
            if (hasTumorRna(case_record) && ! hasTumorDna(case_record)) {
                log.error "targeted mode is not compatible with only tumor RNA provided for ${case_record.case_id}\n\n" +
                    "The targeted workflow requires tumor DNA and can optionally take tumor RNA, depending on " +
                    "the configured panel."
                exit 1
            }

        }

        // Do not allow normal DNA only
        if (hasNormalDna(case_record) && ! hasTumorDna(case_record)) {
            log.error "found only normal DNA input for ${case_record.case_id} but germline only analysis is not supported"
            exit 1
        }

        // Enforce unique samples names within cases
        def sample_ids_duplicated = getSamples(case_record)
            .groupBy { it.sample_id }
            .findAll { k, v -> v.size() > 1 }

        if (sample_ids_duplicated) {
            def duplicate_message_strs = sample_ids_duplicated.collect { sample_id, samples ->
                def key_strs = samples.collect { "${it.sample_type}/${it.sequence_type}" }
                return "  * ${sample_id}: ${key_strs.join(", ")}"
            }
            log.error "duplicate sample names found for ${case_record.case_id}:\n\n${duplicate_message_strs.join("\n")}"
            exit 1
        }

    }

    // NOTE(SW): the following final config checks are performed here since they require additional information
    // regarding processes that are run and also inputs

    def has_alt_contigs = params.genome_type == 'alt'

    // Ensure that custom genomes with ALT contigs that need indexes built have the required .alt file
    def has_bwa_indexes = (params.ref_data_genome_bwamem2_index && params.ref_data_genome_gridss_index)
    def has_alt_file = params.containsKey('ref_data_genome_alt') && params.ref_data_genome_alt
    def run_bwa_or_gridss_index = run_config.stages.alignment && run_config.has_dna_fastq && ! has_bwa_indexes

    if (run_bwa_or_gridss_index && has_alt_contigs && ! has_alt_file) {
        log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  The genome .alt file is required when building bwa-mem2 or GRIDSS indexes\n" +
            "  for reference genomes containing ALT contigs\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

    // Refuse to create STAR index for reference genome containing ALTs, refer to Slack channel
    def run_star_index = run_config.stages.alignment && run_config.has_rna_fastq && ! params.ref_data_genome_star_index

    if (run_star_index && has_alt_contigs) {
        log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Refusing to create the STAR index for a reference genome with ALT contigs.\n" +
            "  Please review https://github.com/alexdobin/STAR docs or contact us on Slack.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

    // Require that an input GTF file is provided when creating STAR index
    if (run_star_index && ! params.ref_data_genome_gtf) {
        log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Creating a STAR index requires the appropriate genome transcript annotations\n" +
            "  as a GTF file. Please contact us on Slack for further information.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

    // Require --isofox_gene_ids argument to be provided in PANEL_RESOURCE_CREATION when RNA inputs are present
    if (run_config.mode == RunMode.PANEL_RESOURCE_CREATION && run_config.has_rna && ! params.isofox_gene_ids) {
        log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Running the panel resource creation workflow with RNA requires that the\n" +
            "  --isofox_gene_ids argument is set with an appropriate input file.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

    // Require panel definition inputs in PANEL_RESOURCE_CREATION mode
    if (run_config.mode == RunMode.PANEL_RESOURCE_CREATION && ! params.target_regions_bed) {
        log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Running the panel resource creation workflow requires that the\n" +
            "  --target_regions_bed argument is set with an appropriate panel target regions file.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

    if (run_config.mode == RunMode.PANEL_RESOURCE_CREATION && ! params.driver_gene_panel) {
        log.error "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Running the panel resource creation workflow requires that the\n" +
            "  --driver_gene_panel argument is set with an appropriate driver gene list file.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        exit 1
    }

}
