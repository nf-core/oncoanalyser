//
// This file holds several functions specific to the workflow/oncoanalyser.nf in the nf-core/oncoanalyser pipeline
//

import static groovy.io.FileType.FILES

import nextflow.Channel
import nextflow.Nextflow

class WorkflowOncoanalyser {

    public static groupByMeta(Map named_args, ... channels) {
        // This method takes multiple channels with the same meta:
        //
        //   channelA: [meta, 'a1']
        //   channelB: [meta, 'b1', 'b2']
        //
        // and combines them into a single channel with the (default) structure:
        //   channel:  [meta, 'a1', 'b1', 'b2']

        // Tag each input channel with its argument position so grouped values can be
        // restored to the original channel order after asynchronous mix().
        def channels_indexed = channels.toList().withIndex().collect { ch, index ->
            ch.map { data ->
                def meta = data[0]
                def values = data[1..-1]
                [meta, [index: index, values: values]]
            }
        }

        // After merging - arrival order after mix() can be random:
        //   channel:
        //     [meta, [index: 1, values: ['b1', 'b2']]]
        //     [meta, [index: 0, values: ['a1']]]
        def channels_merged = Channel.empty().mix(*channels_indexed)

        // After grouping:
        //   channel: [meta, [[index: 1, values: ['b1', 'b2']], [index: 0, values: ['a1']]]]
        def channels_grouped = channels_merged.groupTuple(size: channels.size())

        // After sorting:
        //   channel: [meta, ['a1'], ['b1', 'b2']]]
        def channel_sorted = channels_grouped.map { data ->
            def meta = data[0]
            def values_map = data[1]

            def values_list = values_map
                .sort(false) { it.index }
                .collect { it.values }
            return [meta, *values_list]
        }

        def channel_flat
        def flatten_mode = named_args.getOrDefault('flatten_mode', 'non_recursive')

        // Before: [meta, ['a1'], ['b1', 'b2'], ['c1', ['c2', 'c3']]]
        if (flatten_mode == 'non_recursive') {
            // Flatten one level only
            // After: [meta, 'a1', 'b1', 'b2', 'c1', ['c2', 'c3']]
            channel_flat = channel_sorted.map { data ->
                def meta = data[0]
                def inputs = data[1..-1].collectMany { it }
                return [meta, *inputs]
            }

        } else if (flatten_mode == 'recursive') {
            // Flatten all levels
            // After: [meta, 'a1', 'b1', 'b2', 'c1', 'c2', 'c3']
            channel_flat = channel_sorted.map { it.flatten() }

        } else if(flatten_mode == 'singletons_only') {
            // Flatten top-level singleton lists; preserve longer lists
            // After: [meta, 'a1', ['b1', 'b2'], ['c1', ['c2', 'c3']]]
            channel_flat = channel_sorted.map { data ->
                return data.collect {
                    def is_singleton = it instanceof List && it.size() == 1
                    return is_singleton ? it[0] : it
                }
            }

        } else if(flatten_mode == 'none') {
            channel_flat = channel_sorted

        } else {
            throw new IllegalArgumentException("`flatten_mode` must be one of: 'recursive', 'non_recursive', 'singletons_only', 'none'")
        }

        return channel_flat
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
