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
                    .sort(false) { it.position }
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

    public static getInput(Map named_args, meta, key) {

        def result
        def (key_filetype, key_filetypes, key_sequencetypes) = key

        for (key_sample in [key_filetypes, key_sequencetypes].combinations()) {
            if (meta.containsKey(key_sample) && meta[key_sample].containsKey(key_filetype)) {
                // NOTE(SW): could return early here then false below
                return meta[key_sample].getAt(key_filetype)
                break
            }
        }

        if (result) {
            return result
        } else {
            return false
        }



        /*
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
        */

    }

    // NOTE(SW): function signature required to catch where no named arguments are passed
    public static getInput(ch, key) {
        return getInput([:], ch, key)
    }


    public static joinMeta(Map named_args, ch_a, ch_b) {
        // NOTE(SW): the cross operator is used to allow many-to-one relationship between ch_output
        // and ch_metas
        def key_a = named_args.getOrDefault('key_a', 'group_id')
        def key_b = named_args.getOrDefault('key_b', 'key')

        // Prepare channels
        def ch_ready_a = ch_a.map {
            // We must handle cases where only a process-specific meta is provided, this occurs during process skipping
            // and a meta restore is required to create placeholders in the output.

            if (it instanceof HashMap) {
                return it.getAt(key_b)
            } else {
                return [it[0].getAt(key_b), it[1..-1]]
            }

        }
        def ch_ready_b = ch_b.map { meta -> [meta.getAt(key_a), meta] }

        return ch_ready_b
            .cross(ch_ready_a)
            .map { b, a ->

                // As above, we must accommodate process-specific meta only channel elements, and we would expect just
                // the key to be present for such cases. Otherwise we expect the same key plus additional values in a
                // Collection object.

                if (a instanceof Collection) {
                    def (ka, values) = a
                    def (kb, meta) = b
                    return [meta, *values]
                } else {
                    def (kb, meta) = b
                    return meta
                }

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
