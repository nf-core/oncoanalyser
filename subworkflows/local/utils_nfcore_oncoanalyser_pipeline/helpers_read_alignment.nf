//
// Helpers for read alignment workflows
//

def getFastqsBySampleType(ch_fastq, sample_type) {
    // runnable: channel: [ meta, fastq_info, fastq_fwd, fastq_rev ]
    // skip: channel: [ meta ]
    return ch_fastq
        .branch { meta, fastq_info, fastq_fwd, fastq_rev ->
            def has_inputs = fastq_fwd && fastq_rev
            runnable: fastq_info.sample_type == sample_type && has_inputs
            skip: fastq_info.sample_type == sample_type && ! has_inputs
              return meta
        }
}

def createFastqInputs(ch_inputs_sorted) {
    // Input is a list of the branched channels returned by getFastqsBySampleType
    // Output: channel: [ meta_fastq, fastq_fwd, fastq_rev ]

    def ch_runnables = ch_inputs_sorted.collect { ch_sorted -> ch_sorted.runnable }

    return channel.empty()
        .mix(*ch_runnables)
        .map { meta, fastq_info, fastq_fwd, fastq_rev ->

            // NOTE(SW): initial map sets defaults and conventional ordering of selected fields, merging then overwrites / adds while preserving order
            def rg_id = [
                fastq_info.sample_id,
                fastq_info.library_id,
                fastq_info.lane,
                fastq_info.flowcell
            ].findAll().join('.')

            def rg_entries = [
                ID: rg_id,
                SM: fastq_info.sample_id,
                LB: fastq_info.library_id
            ] + fastq_info.rg_fields

            def rg_line = '@RG\\t' + rg_entries.collect { k, v -> "${k}:${v}" }.join('\\t')

            def meta_fastq = [
                key: meta.group_id,
                id: "${meta.group_id}_${fastq_info.sample_id}_${fastq_info.library_id}_${fastq_info.lane}",
                rg_line: rg_line,
                sample_id: fastq_info.sample_id,
                library_id: fastq_info.library_id,
                lane: fastq_info.lane,
                output_file_id: rg_id,
                sample_type: fastq_info.sample_type,
            ]

            if (fastq_info.flowcell) {
                meta_fastq.id = "${meta_fastq.id}_${fastq_info.flowcell}"
            }

            return [meta_fastq, fastq_fwd, fastq_rev]
        }
}


//
// fastp output channels
//

def expandSplitFastqs(ch_fastp_split) {
    // Input:  channel: [ meta_fastq, [fwd_1, fwd_2, ...], [rev_1, rev_2, ...] ]
    // Output: channel: [ meta_fastq_ready, fastq_fwd, fastq_rev ], one per chunk
    // NOTE(LN): the transpose operator pairs the R1 and R2 chunks by index, and also covers the single chunk case
    // where fastp emits one file per read rather than a list
    return ch_fastp_split
        .transpose()
        .map { meta_fastq, fwd, rev ->

            def split_fwd = fwd.name.replaceAll('\\..+$', '')
            def split_rev = rev.name.replaceAll('\\..+$', '')

            assert split_fwd == split_rev

            // NOTE(SW): split allows meta_fastq_ready to be unique, which is required during reunite below
            def meta_fastq_ready = meta_fastq + [id: "${meta_fastq.id}_${split_fwd}", split: split_fwd]

            return [meta_fastq_ready, fwd, rev]
        }
}

def markUnsplitFastqs(ch_fastq_inputs) {
    // Input:  channel: [ meta_fastq, fastq_fwd, fastq_rev ]
    // Output: channel: [ meta_fastq_ready, fastq_fwd, fastq_rev ]
    return ch_fastq_inputs
        .map { meta_fastq, fastq_fwd, fastq_rev ->
            def meta_fastq_ready = meta_fastq + [split: null]
            return [meta_fastq_ready, fastq_fwd, fastq_rev]
        }
}


//
// BWA-MEM2 channels
//

def createBwamem2Inputs(ch_fastqs_ready) {
    // Input:  channel: [ meta_fastq_ready, fastq_fwd, fastq_rev ]
    // Output: channel: [ meta_bwamem2, fastq_fwd, fastq_rev ]
    return ch_fastqs_ready
        .map { meta_fastq_ready, fastq_fwd, fastq_rev ->
            def meta_bwamem2 = meta_fastq_ready.clone()
            return [meta_bwamem2, fastq_fwd, fastq_rev]
        }
}

def getSampleFastqCounts(ch_bwamem2_inputs) {
    // Count expected BAMs per sample for non-blocking groupTuple op
    // Input:  channel: [ meta_bwamem2, fastq_fwd, fastq_rev ]
    // Output: channel: [ meta_group, group_size ]
    return ch_bwamem2_inputs
        .map { meta_bwamem2, _reads_fwd, _reads_rev ->

            def meta_group = [
                key: meta_bwamem2.key,
                sample_type: meta_bwamem2.sample_type,
            ]

            return [meta_group, meta_bwamem2]
        }
        .groupTuple()
        .map { meta_group, metas_bwamem2 -> return [meta_group, metas_bwamem2.size()] }
}
