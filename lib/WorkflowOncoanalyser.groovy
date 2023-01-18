//
// This file holds several functions specific to the workflow/oncoanalyser.nf in the nf-core/oncoanalyser pipeline
//

import static groovy.io.FileType.FILES

import nextflow.Channel
import nextflow.Nextflow

import Constants

class WorkflowOncoanalyser {

    //
    // Check and validate parameters
    //
    public static void initialise(params, workflow, log) {

        // TODO(SW): allow users to set all appropriate reference genomes manually in config or CLI, including version and type (see below)

        if (!params.genome) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome must be set using the --genome CLI argument or in a configuration file.\n" +
                "  Currently, the available genome are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        } else if (!params.genomes.containsKey(params.genome)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        // NOTE(SW): restricting allowable genome values to GRCh37_hmf for now
        if (params.genome != 'GRCh37_hmf') {
            log.error "ERROR: currently only the GRCh37_hmf genome is supported but got \"${params.genome}\"" +
                ", please adjust the --genome argument accordingly."
            System.exit(1)
        }

        if (Constants.GENOMES_VERSION_37.contains(params.genome)) {
            params.ref_data_genome_version = '37'
        } else if (Constants.GENOMES_VERSION_38.contains(params.genome)) {
            params.ref_data_genome_version = '38'
        } else {
            def genome_version_list_all = Constants.GENOMES_VERSION_37 + Constants.GENOMES_VERSION_38
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' is not defined in genome version list.                 \n" +
                "  Currently, the list of genomes in the version list include:\n" +
                "  ${genome_version_list_all.join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }

        // TODO(SW): when allowing user to set custom genome, require this to be explicitly set
        if (Constants.GENOMES_ALT.contains(params.genome)) {
            params.ref_data_genome_type = 'alt'
        } else {
            params.ref_data_genome_type = 'no_alt'
        }

        if (!params.containsKey('ref_data_hmf_data_base')) {
            if (params.ref_data_genome_version == '37') {
                params.ref_data_hmf_data_base = Constants.HMF_DATA_37_BASE
            } else if (params.ref_data_genome_version == '38') {
                params.ref_data_hmf_data_base = Constants.HMF_DATA_38_BASE
            } else {
                assert false : "Got a bad genome version: ${params.ref_data_genome_version}"
            }
        }

        if (!params.containsKey('ref_data_virusbreakenddb_path')) {
            params.ref_data_virusbreakenddb_path = Constants.VIRUSBREAKENDDB_PATH
        }

        def null_check = [
           'ref_data_genome_fasta',
           'ref_data_genome_type',
           'ref_data_genome_version',
           'ref_data_virusbreakenddb_path',
        ]
        null_check.each { k ->
            if (!params[k]) {
                log.error "ERROR: '${k}' cannot be set to null in any configuration and must be adjusted or removed to proceed."
                System.exit(1)
            }
        }

        // Download region file for collectwgsmetrics if in test mode
        // NOTE(SW): this will be removed as part of the overhaul for testing
        if (workflow.profile.contains('test')) {
            def stage_dir = new File(workflow.workDir.toString(), 'stage/manual/')
            def interval_file = new File(stage_dir, 'collectwgsmetrics.interval_list')
            if (! interval_file.exists()) {
                stage_dir.mkdirs()
                interval_file.createNewFile()
                interval_file << new URL (params.ref_data_wgsmetrics_intervals_url).getText()
            }
            params.ref_data_wgsmetrics_intervals_local = interval_file
        }
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

        if (named_args.get('flatten', true)) {
            def flatten_mode = named_args.get('flatten_mode', 'recursive')
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

    public static getInput(ch, key) {
        return ch.map { meta -> [meta, meta.getAt(key)] }
    }

    public static joinMeta(Map named_args, ch_a, ch_b) {
        // NOTE(SW): the cross operator is used to allow many-to-one relationship between ch_output
        // and ch_metas
        def key_a = named_args.get('key_a', 'id')
        def key_b = named_args.get('key_b', 'key')
        def ch_ready_a = ch_a.map { [it[0].get(key_b), it[1..-1]] }
        def ch_ready_b = ch_b.map { meta -> [meta.get(key_a), meta] }
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
