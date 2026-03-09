//
// This file holds several functions specific to the workflow/oncoanalyser.nf in the nf-core/oncoanalyser pipeline
//

import static groovy.io.FileType.FILES

import nextflow.Channel
import nextflow.Nextflow

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

        def flatten_mode = named_args.getOrDefault('flatten_mode', 'non_recursive')
        if (flatten_mode == 'recursive') {
            return r.map { it.flatten() }
        } else if (flatten_mode == 'non_recursive') {
            return r.map { data ->
                def meta = data[0]
                def inputs = data[1..-1].collectMany { it }
                return [meta, *inputs]
            }
        } else if(flatten_mode == 'singletons_only') {
            return r.map { data ->
                return data.collect {
                    def is_singleton = it instanceof List && it.size() == 1
                    return is_singleton ? it[0] : it
                }
            }
        } else if(flatten_mode == 'none') {
            return r
        } else {
            throw new IllegalArgumentException("`flatten_mode` must be one of: 'recursive', 'non_recursive', 'singletons_only', 'none'")
        }
    }

    // NOTE(SW): function signature required to catch where no named arguments are passed
    public static groupByMeta(... channels) {
        return groupByMeta([:], *channels)
    }

    public static restoreMeta(ch_output, ch_metas) {
        // NOTE(SW): ch_output must contain a Map in the first position with a key named 'key' that
        // contains the corresponding meta.id value, for example: [val(meta_process), *process_outputs]

        def ch_ready_a = ch_output.map { [it[0].getAt('key'), it[1..-1]] }
        def ch_ready_b = ch_metas.map { meta -> [meta.getAt('group_id'), meta] }
        return ch_ready_b
            // NOTE(SW): the cross operator is used to allow many-to-one relationship between ch_output
            // and ch_metas
            .cross(ch_ready_a)
            .map { b, a ->
                def (ka, values) = a
                def (kb, meta) = b
                return [meta, *values]
            }
    }
}
