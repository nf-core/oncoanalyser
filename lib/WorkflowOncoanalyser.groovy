//
// This file holds several functions specific to the workflow/oncoanalyser.nf in the nf-core/oncoanalyser pipeline
//

import static groovy.io.FileType.FILES

import nextflow.Channel
import nextflow.Nextflow

import Constants
import Processes
import Utils

class WorkflowOncoanalyser {

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
        if (params.containsKey('run_mode')) {
            run_mode = Utils.getRunMode(params.run_mode, log)
        } else {
            // Bad configuration, catch in validateParams
            return
        }

        if (run_mode == Constants.RunMode.WTS) {
            if (!params.containsKey('run_type')) {
                params.run_type = 'tumor_only'
            }
        }

        def stages = Processes.getRunStages(
            run_mode,
            params.processes_include,
            params.processes_exclude,
            params.processes_manual,
            log,
        )

        if (stages.virusinterpreter) {
            if (!params.containsKey('ref_data_virusbreakenddb_path')) {
                params.ref_data_virusbreakenddb_path = Constants.VIRUSBREAKENDDB_PATH
            }
        }

        if (stages.lilac) {
            if (params.ref_data_genome_version == '38' && params.ref_data_genome_type == 'alt' && !params.containsKey('ref_data_hla_slice_bed')) {
                params.ref_data_hla_slice_bed = Constants.HLA_SLICE_BED_GRCH38_ALT_PATH
            }
        }
    }

    //
    // Check and validate parameters
    //
    public static void validateParams(params, log) {

        // Common parameters

        if (!params.genome) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
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
                "  genome version list.\n" +
                "  Currently, the list of genomes in the version list include:\n" +
                "  ${Constants.GENOMES_DEINFED.join(", ")}\n" +
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

        if (!params.run_mode) {
            def run_modes = Utils.getEnumNames(Constants.RunMode).join('\n    - ')
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Run mode must be set using the --run_mode CLI argument or in a configuration  \n" +
                "  file.\n" +
                "  Currently, the available run modes are:\n" +
                "    - ${run_modes}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        if (!params.run_type) {
            def run_types = Utils.getEnumNames(Constants.RunType).join('\n    - ')
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Run type must be set using the --run_type CLI argument or in a configuration  \n" +
                "  file.\n" +
                "  Currently, the available run types are:\n" +
                "    - ${run_types}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        def run_mode = Utils.getRunMode(params.run_mode, log)
        def run_type = Utils.getRunType(params.run_type, log)

        if (run_mode == Constants.RunMode.WTS) {
            if (run_type != Constants.RunType.TUMOR_ONLY) {
                log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  The WTS run mode does not support tumor/normal data, please adjust the CLI \n" +
                    "  --run_type option or corresponding configuration file.\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                System.exit(1)
            }

        }

    }

    public static String paramsSummaryLog(workflow, params, log) {
        def summary_log = ''
        summary_log += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        summary_log += NfcoreSchema.paramsSummaryLog(workflow, params)
        summary_log += '\n' + WorkflowMain.citation(workflow) + '\n'
        summary_log += NfcoreTemplate.dashedLine(params.monochrome_logs)
        log.info summary_log
    }

    public static getRunConfig(params, log) {
        def run_mode = Utils.getRunMode(params.run_mode, log)
        def run_type = Utils.getRunType(params.run_type, log)
        def stages = Processes.getRunStages(
            run_mode,
            params.processes_include,
            params.processes_exclude,
            params.processes_manual,
            log,
        )

        return [
            mode: run_mode,
            type: run_type,
            stages: stages,
        ]
    }

    public static groupByMeta(Map named_args, ... channels) {
        def r = channels
        // Set position; required to use non-blocking .mix operator
        // NOTE(SW): operating on native list object containing channels
        def i = 0
        r = r
            .collect { ch ->
                def ii = i
                def d = ch.map { data ->
                    def meta = data[0]
                    def values = data[1..-1]
                    return [meta, [position: ii, values: values]]
                }
                i++
                return d
            }

        r = Channel.empty().mix(*r)

        // NOTE(SW): As of Nextflow 22.10.6, groupTuple requires a matching meta /and/ an additional element to complete without error, these placeholders are filtered in the groupByMeta function
        r = r.filter { it[0] != Constants.PLACEHOLDER_META }

        r = r
            .groupTuple(size: channels.size())
            .map { data ->
                def meta = data[0]
                def values_map = data[1]

                def values_list = values_map
                    .sort { it.position }
                    .collect { it.values }
                return [meta, *values_list]
            }

        if (named_args.getOrDefault('flatten', true)) {
            def flatten_mode = named_args.getOrDefault('flatten_mode', 'nonrecursive')
            if (flatten_mode == 'recursive') {
                r = r.map { it.flatten() }
            } else if (flatten_mode == 'nonrecursive') {
                r = r.map { data ->
                    def meta = data[0]
                    def inputs = data[1..-1].collectMany { it }
                    return [meta, *inputs]
                }
            } else {
                System.err.println "ERROR: got bad flatten_mode: ${flatten_mode}"
                System.exit(1)
            }
        }

        return r
    }

    // NOTE(SW): function signature required to catch where no named arguments are passed
    public static groupByMeta(... channels) {
        return groupByMeta([:], *channels)
    }

    public static getInput(Map named_args, ch, key) {
        def input_type = named_args.getOrDefault('type', 'required')
        return ch.map { meta ->
            def result
            for (key_combination in key.combinations()) {
                if (meta.containsKey(key_combination)) {
                    result = [meta, meta.getAt(key_combination)]
                    break
                }
            }

            if (result) {
                return result
            }

            if (input_type == 'required') {
                return [Constants.PLACEHOLDER_META, null]
            } else if (input_type == 'optional') {
                return [meta, []]
            } else {
                System.err.println "ERROR: got bad input type: ${input_type}"
                System.exit(1)
            }
        }
    }

    // NOTE(SW): function signature required to catch where no named arguments are passed
    public static getInput(ch, key) {
        return getInput([:], ch, key)
    }


    public static joinMeta(Map named_args, ch_a, ch_b) {
        // NOTE(SW): the cross operator is used to allow many-to-one relationship between ch_output
        // and ch_metas
        def key_a = named_args.getOrDefault('key_a', 'id')
        def key_b = named_args.getOrDefault('key_b', 'key')
        def ch_ready_a = ch_a.map { [it[0].getAt(key_b), it[1..-1]] }
        def ch_ready_b = ch_b.map { meta -> [meta.getAt(key_a), meta] }
        return ch_ready_b
            .cross(ch_ready_a)
            .map { b, a ->
                def (ka, values) = a
                def (kb, meta) = b
                return [meta, *values]
            }
    }

    // NOTE(SW): function signature required to catch where no named arguments are passed
    public static joinMeta(ch_output, ch_metas) {
        joinMeta([:], ch_output, ch_metas)
    }

    public static restoreMeta(ch_output, ch_metas) {
        // NOTE(SW): ch_output must contain a Map in the first position with a key named 'key' that
        // contains the corresponding meta.id value, for example: [val(meta_process), *process_outputs]
        joinMeta([:], ch_output, ch_metas)
    }
}
