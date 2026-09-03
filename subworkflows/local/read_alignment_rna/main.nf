//
// Align RNA reads
//

include { BWAMEM2_ALIGN_RNA } from '../../../modules/local/bwa-mem2/mem/rna/main'
include { FASTP_SPLIT       } from '../../../modules/local/fastp/split/main'

include { getFastqsBySampleType } from '../utils_nfcore_oncoanalyser_pipeline/helpers_read_alignment'
include { createFastqInputs     } from '../utils_nfcore_oncoanalyser_pipeline/helpers_read_alignment'
include { expandSplitFastqs     } from '../utils_nfcore_oncoanalyser_pipeline/helpers_read_alignment'
include { markUnsplitFastqs     } from '../utils_nfcore_oncoanalyser_pipeline/helpers_read_alignment'
include { createBwamem2Inputs   } from '../utils_nfcore_oncoanalyser_pipeline/helpers_read_alignment'
include { getSampleFastqCounts  } from '../utils_nfcore_oncoanalyser_pipeline/helpers_read_alignment'

workflow READ_ALIGNMENT_RNA {
    take:
    // Sample data
    ch_inputs            // channel: [mandatory] [ meta ]
    ch_fastq             // channel: [mandatory] [ meta, fastq_info, fastq_fwd, fastq_rev ]

    // Reference data
    genome_fasta         // channel: [mandatory] /path/to/genome_fasta
    genome_bwamem2_index // channel: [mandatory] /path/to/genome_bwa-mem2_index_dir/

    // Params
    max_fastq_records    // numeric: [optional]  max number of FASTQ records per split

    main:
    //
    // STEP: Handle inputs
    //
    // Sort inputs
    // runnable: channel: [ meta, fastq_info, fastq_fwd, fastq_rev ]
    // skip: channel: [ meta ]

    // NOTE(LN): RNA is only ever sequenced for the tumor sample
    ch_inputs_rna_sorted = getFastqsBySampleType(ch_fastq, 'tumor')

    // Create FASTQ input channel
    // channel: [ meta_fastq, fastq_fwd, fastq_rev ]
    ch_fastq_inputs = createFastqInputs([ch_inputs_rna_sorted])

    //
    // MODULE: fastp
    //
    // Split FASTQ into chunks if requested for distributed processing
    // channel: [ meta_fastq_ready, fastq_fwd, fastq_fwd ]
    ch_fastqs_ready = channel.empty()
    // NOTE(SW): required for strict syntax without params block declaration
    if (max_fastq_records.toInteger() > 0) {

        // Run process
        FASTP_SPLIT(
            ch_fastq_inputs,
            // NOTE(SW): required for strict syntax without params block declaration
            max_fastq_records.toInteger(),
        )

        ch_fastqs_ready = expandSplitFastqs(FASTP_SPLIT.out[0])

    } else {

        ch_fastqs_ready = markUnsplitFastqs(ch_fastq_inputs)

    }

    //
    // MODULE: BWA-MEM2
    //
    // Create process input channel
    // channel: [ meta_bwamem2, fastq_fwd, fastq_rev ]
    ch_bwamem2_inputs = createBwamem2Inputs(ch_fastqs_ready)

    // Run process
    BWAMEM2_ALIGN_RNA(
        ch_bwamem2_inputs,
        genome_fasta,
        genome_bwamem2_index,
    )

    // Reunite BAMs
    // channel: [ meta_group, group_size ]
    ch_sample_fastq_counts = getSampleFastqCounts(ch_bwamem2_inputs)

    // Now, group with expected size
    // NOTE(LN): RNA alignments are name-grouped and therefore unindexed, so no index is carried here
    // channel: [ meta_group, [aln, ...] ]
    ch_alns_united = ch_sample_fastq_counts
        // channel: [ [ meta_group, count ], [ meta_group, aln ] ]
        .cross(
            // First element to match meta_group above for `cross`
            channel.topic('bwamem2_align_rna_bam').map { meta_bwamem2, aln -> [[key: meta_bwamem2.key, sample_type: meta_bwamem2.sample_type], aln] }
        )
        .map { count_tuple, inputs_tuple ->
            def group_size = count_tuple[1]
            def (meta_group, aln) = inputs_tuple

            return tuple(groupKey(meta_group, group_size), aln)
        }
        .groupTuple()

    //
    // STEP: Handle outputs
    //
    // Set outputs, restoring original meta
    // NOTE(LN): the empty index position is kept so that the output shape matches the DNA alignment subworkflow and
    // the placeholders set by the calling workflows
    // channel: [ meta, [aln, ...], [] ]
    ch_outputs_rna = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_alns_united.map { meta_group, alns -> [meta_group, alns, []] }, ch_inputs),
            ch_inputs_rna_sorted.skip.unique().map { meta -> [meta, [], []] },
        )

    emit:
    rna = ch_outputs_rna // channel: [ meta, [aln, ...], [] ]
}
