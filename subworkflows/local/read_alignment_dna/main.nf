//
// Align DNA reads
//

include { BWAMEM2_ALIGN_DNA } from '../../../modules/local/bwa-mem2/mem/dna/main'
include { FASTP_SPLIT       } from '../../../modules/local/fastp/split/main'

include { getFastqsBySampleType } from '../utils_nfcore_oncoanalyser_pipeline/helpers_read_alignment'
include { createFastqInputs     } from '../utils_nfcore_oncoanalyser_pipeline/helpers_read_alignment'
include { expandSplitFastqs     } from '../utils_nfcore_oncoanalyser_pipeline/helpers_read_alignment'
include { markUnsplitFastqs     } from '../utils_nfcore_oncoanalyser_pipeline/helpers_read_alignment'
include { createBwamem2Inputs   } from '../utils_nfcore_oncoanalyser_pipeline/helpers_read_alignment'
include { getSampleFastqCounts  } from '../utils_nfcore_oncoanalyser_pipeline/helpers_read_alignment'

workflow READ_ALIGNMENT_DNA {
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

    // NOTE(SW): the logic block diverges from the standard to accommodate the nature of inputs
    //   - existing input inferred from upstream check
    //   - channel elements are each a single FASTQ pair
    //   - hence case meta maps to pairs as one-to-many
    //   - skip conditional is adjusted for this albiet emits duplicate case meta
    //   - case meta in skip channel deduplicated at join below

    ch_inputs_tumor_sorted = getFastqsBySampleType(ch_fastq, 'tumor')
    ch_inputs_normal_sorted = getFastqsBySampleType(ch_fastq, 'normal')
    ch_inputs_donor_sorted = getFastqsBySampleType(ch_fastq, 'donor')

    // Create FASTQ input channel
    // channel: [ meta_fastq, fastq_fwd, fastq_rev ]
    ch_fastq_inputs = createFastqInputs([
        ch_inputs_tumor_sorted,
        ch_inputs_normal_sorted,
        ch_inputs_donor_sorted,
    ])

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
    BWAMEM2_ALIGN_DNA(
        ch_bwamem2_inputs,
        genome_fasta,
        genome_bwamem2_index,
    )

    // Reunite BAMs
    // channel: [ meta_group, group_size ]
    ch_sample_fastq_counts = getSampleFastqCounts(ch_bwamem2_inputs)

    // Now, group with expected size then sort into tumor and normal channels
    // channel: [ meta_group, [aln, ...], [idx, ...] ]
    ch_alns_united = ch_sample_fastq_counts
        // channel: [ [ meta_group, count ], [ meta_group, aln, idx ] ]
        .cross(
            // First element to match meta_group above for `cross`
            channel.topic('bwamem2_align_dna_bam').map { meta_bwamem2, aln, idx -> [[key: meta_bwamem2.key, sample_type: meta_bwamem2.sample_type], aln, idx] }
        )
        .map { count_tuple, inputs_tuple ->
            def group_size = count_tuple[1]
            def (meta_group, aln, idx) = inputs_tuple

            return tuple(groupKey(meta_group, group_size), aln, idx)
        }
        .groupTuple()
        .branch { meta_group, alns, idxs ->
            assert ['tumor', 'normal', 'donor'].contains(meta_group.sample_type)
            tumor: meta_group.sample_type == 'tumor'
            normal: meta_group.sample_type == 'normal'
            donor: meta_group.sample_type == 'donor'
            placeholder: true
        }

    //
    // STEP: Handle outputs
    //
    // Set outputs, restoring original meta
    // channel: [ meta, [aln, ...], [idx, ...] ]
    ch_outputs_tumor = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_alns_united.tumor, ch_inputs),
            ch_inputs_tumor_sorted.skip.unique().map { meta -> [meta, [], []] },
        )

    // channel: [ meta, [aln, ...], [idx, ...] ]
    ch_outputs_normal = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_alns_united.normal, ch_inputs),
            ch_inputs_normal_sorted.skip.unique().map { meta -> [meta, [], []] },
        )

    // channel: [ meta, [aln, ...], [idx, ...] ]
    ch_outputs_donor = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_alns_united.donor, ch_inputs),
            ch_inputs_donor_sorted.skip.unique().map { meta -> [meta, [], []] },
        )

    emit:
    tumor  = ch_outputs_tumor  // channel: [ meta, [aln, ...], [idx, ...] ]
    normal = ch_outputs_normal // channel: [ meta, [aln, ...], [idx, ...] ]
    donor  = ch_outputs_donor  // channel: [ meta, [aln, ...], [idx, ...] ]
}
