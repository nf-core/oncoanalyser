//
// Channel helpers for the nf-core/oncoanalyser pipeline
//

include { getDonorDnaSamples  } from './accessors_samples'
include { getNormalDnaSample  } from './accessors_samples'
include { getTumorDnaSample   } from './accessors_samples'
include { getTumorRnaSample   } from './accessors_samples'
include { hasDonorDnaAlns     } from './accessors_alignments'
include { hasDonorDnaFastqs   } from './accessors_samples'
include { hasInput            } from './accessors_samples'
include { hasNormalDnaAln     } from './accessors_alignments'
include { hasNormalDnaFastq   } from './accessors_samples'
include { hasTumorDnaAln      } from './accessors_alignments'
include { hasTumorDnaFastq    } from './accessors_samples'
include { hasTumorRnaFastq    } from './accessors_samples'
include { FileType            } from './types_enums'

def getPlaceholderMeta() { return [meta_placeholder: null] }

def groupByMeta(Map named_args, List channels) {
    def r = channels
    // Set position; required to use non-blocking .mix operator
    // NOTE(SW): operating on native list object containing channels
    def i = 0
    r = r.collect { ch ->
        def ii = i
        def d = ch.map { data ->
            def meta = data[0]
            def values = data[1..-1]
            return [meta, [position: ii, values: values]]
        }
        i = i + 1
        return d
    }

    // NOTE(SW): seed the fold with the first channel, not Channel.empty(). A typed process
    // `topic:` output is a DataflowStream; mixing it with a queue (Channel.empty()) fails in
    // typed mode. Use `.drop(1)` rather than `r[1..-1]`: the reverse range 1..-1 normalises
    // to subList(0, 2) for a single-element list and throws IndexOutOfBoundsException.
    r = r.drop(1).inject(r[0]) { acc, ch -> acc.mix(ch) }

    // NOTE(SW): As of Nextflow 22.10.6, groupTuple requires a matching meta /and/ an additional element to complete without error, these placeholders are filtered in the groupByMeta function
    r = r.filter { it[0] != getPlaceholderMeta() }

    r = r
        .groupTuple(size: channels.size())
        .map { data ->
            def meta = data[0]
            def values_map = data[1]

            def values_list = values_map.sort(false) { it.position }.collect { it.values }
            return [meta] + values_list
        }

    if (named_args.getOrDefault('flatten', true)) {
        def flatten_mode = named_args.getOrDefault('flatten_mode', 'nonrecursive')
        if (flatten_mode == 'recursive') {
            r = r.map { it.flatten() }
        } else if (flatten_mode == 'nonrecursive') {
            r = r.map { data ->
                def meta = data[0]
                def inputs = data[1..-1].collectMany { it }
                return [meta] + inputs
            }
        } else {
            error("got bad flatten_mode: ${flatten_mode}")
        }
    }

    return r
}

// NOTE(SW): function signature required to catch where no named arguments are passed
def groupByMeta(List channels) {
    return groupByMeta([:], channels)
}

def joinMeta(Map named_args, ch_a, ch_b) {
    // NOTE(SW): the cross operator is used to allow many-to-one relationship between ch_output
    // and ch_metas
    def key_a = named_args.getOrDefault('key_a', 'case_id')
    def key_b = named_args.getOrDefault('key_b', 'key')
    def ch_ready_a = ch_a.map { [it[0].getAt(key_b), it[1..-1]] }
    def ch_ready_b = ch_b.map { meta -> [meta.getAt(key_a), meta] }
    return ch_ready_b
        .cross(ch_ready_a)
        .map { b, a ->
            def (ka, values) = a
            def (kb, meta) = b
            return [meta] + values
        }
}

// NOTE(SW): function signature required to catch where no named arguments are passed
def joinMeta(ch_output, ch_metas) {
    joinMeta([:], ch_output, ch_metas)
}

def restoreMeta(ch_output, ch_metas) {
    // NOTE(SW): ch_output must contain a Map in the first position with a key named 'key' that
    // contains the corresponding meta.id value, for example: [val(meta_process), *process_outputs]
    joinMeta([:], ch_output, ch_metas)
}

def getDnaFastqChannel(ch_inputs) {
    // Sort inputs
    // channel: [ case_record ]
    def ch_inputs_tumor_sorted = ch_inputs
        .branch { case_record ->
            def has_existing = hasTumorDnaAln(case_record)
            runnable: hasTumorDnaFastq(case_record) && ! has_existing
            skip: true
        }

    def ch_inputs_normal_sorted = ch_inputs
        .branch { case_record ->
            def has_existing = hasNormalDnaAln(case_record)
            runnable: hasNormalDnaFastq(case_record) && ! has_existing
            skip: true
        }

    def ch_inputs_donor_sorted = ch_inputs
        .branch { case_record ->
            def has_existing = hasDonorDnaAlns(case_record)
            runnable: hasDonorDnaFastqs(case_record) && ! has_existing
            skip: true
        }

    // Create FASTQ input channel
    // channel: [ case_record, fastq_info, fastq_fwd, fastq_rev ]
    def ch_fastqs = channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable.map { case_record -> [case_record, getTumorDnaSample(case_record), 'tumor'] },
            ch_inputs_normal_sorted.runnable.map { case_record -> [case_record, getNormalDnaSample(case_record), 'normal'] },
            ch_inputs_donor_sorted.runnable.flatMap { case_record -> getDonorDnaSamples(case_record).collect { donor -> [case_record, donor, 'donor'] } },
        )
        .flatMap { case_record, sample, sample_type ->
            sample.files.getAt(FileType.FASTQ)
                .collect { fastq ->
                    def fastq_info = [
                        'sample_id': sample.sample_id,
                        'library_id': fastq.library_id,
                        'lane': fastq.lane,
                        'sample_type': sample_type,
                        'single_end': fastq.single_end,
                        'rg_fields': fastq.rg_fields,
                    ]

                    if (fastq.flowcell) {
                         fastq_info.flowcell = fastq.flowcell
                    }

                    return [case_record, fastq_info, fastq.read_fwd, fastq.read_rev ?: null]
                }
        }

    return channel.empty()
        .mix(
            ch_fastqs,
            ch_inputs_tumor_sorted.skip.map { case_record -> [case_record, [:], null, null] },
            ch_inputs_normal_sorted.skip.map { case_record -> [case_record, [:], null, null] },
            ch_inputs_donor_sorted.skip.map { case_record -> [case_record, [:], null, null] },
        )
}

def getRnaFastqChannel(ch_inputs) {
    // Sort inputs
    // channel: [ case_record ]
    def ch_inputs_sorted = ch_inputs
        .branch { case_record ->
            def has_existing = hasInput(getTumorRnaSample(case_record), FileType.ALN)
            runnable: hasTumorRnaFastq(case_record) && ! has_existing
            skip: true
        }

    // Create FASTQ input channel
    // channel: [ case_record, fastq_info, fastq_fwd, fastq_rev ]
    def ch_fastqs = ch_inputs_sorted.runnable
        .flatMap { case_record ->
            def sample = getTumorRnaSample(case_record)
            sample.files
                .getAt(FileType.FASTQ)
                .collect { fastq ->
                    def fastq_info = [
                        'sample_id': sample.sample_id,
                        'library_id': fastq.library_id,
                        'lane': fastq.lane,
                        'rg_fields': fastq.rg_fields,
                    ]

                    if (fastq.flowcell) {
                         fastq_info.flowcell = fastq.flowcell
                    }

                    return [case_record, fastq_info, fastq.read_fwd, fastq.read_rev]
                }
        }

    return channel.empty()
        .mix(
            ch_fastqs,
            ch_inputs_sorted.skip.map { case_record -> [case_record, [:], null, null] },
        )
}
