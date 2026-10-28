//
// Process read UMIs
//

include { FASTP_UMI   } from '../../../modules/local/fastp/umi/main'
include { FASTQ_TOOLS } from '../../../modules/local/fastqtools/main'

workflow READ_UMI_PROCESSING {
    take:
    // Sample data
    ch_inputs               // channel: [mandatory] [ meta ]
    ch_dna_fastq            // channel: [mandatory] [ meta, fastq_info, fastq_fwd, fastq_rev ]
    ch_rna_fastq            // channel: [mandatory] [ meta, fastq_info, fastq_fwd, fastq_rev ]

    // Reference data
    known_umis              // channel: [mandatory] /path/to/known_umis_file

    // Params
    fastp_umi_enabled       // boolean: [mandatory] enable fastp UMI processing
    fastp_umi_location      //  string: [optional]  fastp UMI location argument (--umi_loc)
    fastp_umi_length        // numeric: [optional]  fastp UMI length argument (--umi_len)
    fastp_umi_skip          // numeric: [optional]  fastp UMI skip argument (--umi_skip)
    fastq_tools_umi_enabled // boolean: [mandatory] enable fastq-tools UMI processing
    fastq_tools_umi_delim   // boolean: [optional]  fastq-tools -umi_delim argument

    main:
    //
    // STEP: Handle inputs
    //
    // Sort inputs
    // runnable: channel: [ meta, sequence_type, fastq_info, fastq_fwd, fastq_rev ]
    // skip: channel: [ meta ]
    ch_inputs_dna_sorted = ch_dna_fastq
        .branch { meta, fastq_info, fastq_fwd, fastq_rev ->
            // NOTE(SW): inferred state from upstream
            def has_inputs = fastq_fwd && fastq_rev
            runnable: has_inputs
                return [meta, 'dna', fastq_info, fastq_fwd, fastq_rev]
            skip: true
              return meta
        }

    ch_inputs_rna_sorted = ch_rna_fastq
        .branch { meta, fastq_info, fastq_fwd, fastq_rev ->
            // NOTE(SW): inferred state from upstream
            def has_inputs = fastq_fwd && fastq_rev
            runnable: has_inputs
                return [meta, 'rna', fastq_info, fastq_fwd, fastq_rev]
            skip: true
                return meta
        }

    // Create base FASTQ input channel
    // channel: [ meta, sequence_type, fastq_info, fastq_fwd, fastq_rev ]
    ch_inputs_runnable = channel.empty()
        .mix(
            ch_inputs_dna_sorted.runnable,
            ch_inputs_rna_sorted.runnable,
        )

    // channel: [ meta_fastq, fastq_fwd, fastq_rev ]
    ch_fastq_inputs = ch_inputs_runnable
        .map { meta, sequence_type, fastq_info, fastq_fwd, fastq_rev ->

              def meta_fastq = [
                  key: meta.group_id,
                  id: "${meta.group_id}_${fastq_info.sample_id}",
                  sequence_type: sequence_type,
                  sample_id: fastq_info.sample_id,
                  library_id: fastq_info.library_id,
                  lane: fastq_info.lane,
                  flowcell: fastq_info.flowcell,
                  rg_fields: fastq_info.rg_fields,
              ]

              if (sequence_type == 'dna') {
                  meta_fastq.sample_type = fastq_info.sample_type
              }

              return [meta_fastq, fastq_fwd, fastq_rev]

        }

    // Process UMIs
    // The run conditions for each stage is as follows:
    //  - DNA: either fastp or fastqtools
    //  - RNA: fastqtools only
    //
    // As such DNA / RNA may trigger both fastp (DNA) and fastqtools (RNA), so each must be handled separately

    //
    // MODULE: fastp
    //
    // channel: [ meta_fastq, fastq_fwd, fastq_rev ]
    ch_post_fastp = channel.empty()
    if (fastp_umi_enabled) {


        // Sort inputs
        // channel: runnable: [ meta_fastq, fastq_fwd, fastq_rev ]
        // channel: skip: [ meta_fastq, fastq_fwd, fastq_rev ]
        ch_fastp_inputs_sorted = ch_fastq_inputs
            .branch { meta_fastq, fastq_fwd, fastq_rev ->
                runnable: meta_fastq.sequence_type == 'dna'
                skip: true
            }

        // Run process
        FASTP_UMI(
            ch_fastp_inputs_sorted.runnable,
            fastp_umi_location,
            fastp_umi_length,
            // NOTE(SW): required for strict syntax without params block declaration
            fastp_umi_skip.toInteger(),
        )


        // Set outputs
        ch_post_fastp = channel.empty()
            .mix(
                channel.topic('fastp_umi_fastq'),
                ch_fastp_inputs_sorted.skip,
            )

    } else {

        ch_post_fastp = ch_fastq_inputs

    }

    //
    // MODULE: FASTQTOOLS
    //
    // channel: [ meta_fastq, fastq_fwd, fastq_rev ]
    ch_post_fastqtools = channel.empty()
    if (fastq_tools_umi_enabled) {

        // NOTE(SW): only run DNA FASTQs when fastp hasn't already been run
        // Sort inputs
        // channel: runnable: [ meta_fastq, fastq_fwd, fastq_rev ]
        // channel: skip: [ meta_fastq, fastq_fwd, fastq_rev ]
        ch_fastqtools_inputs_sorted = ch_post_fastp
            .branch { meta_fastq, fastq_fwd, fastq_rev ->
                runnable: ! (fastp_umi_enabled && meta_fastq.sequence_type == 'dna')
                skip: true
            }

        // Run process
        FASTQ_TOOLS(
            ch_fastqtools_inputs_sorted.runnable,
            fastq_tools_umi_delim,
            known_umis,
        )

        // Set outputs
        ch_post_fastqtools = channel.empty()
            .mix(
                channel.topic('fastqtools_fastq'),
                ch_fastqtools_inputs_sorted.skip,
            )

    } else {

        ch_post_fastqtools = ch_post_fastp

    }

    // Re-construct fastq_info and separate processed FASTQ into DNA / RNA sequence type
    // NOTE(SW): not taking the route of grouping since identity requires additional keying in this one-to-many scenario
    ch_fastq_processed_sorted = ch_post_fastqtools
        .map { meta_fastq, fastq_fwd, fastq_rev ->

            def fastq_info = [
                'sample_id': meta_fastq.sample_id,
                'library_id': meta_fastq.library_id,
                'lane': meta_fastq.lane,
                'flowcell': meta_fastq.flowcell,
                'rg_fields': meta_fastq.rg_fields,
            ]

            if (meta_fastq.sequence_type == 'dna') {
                fastq_info.sample_type = meta_fastq.sample_type
            }

            return [meta_fastq, fastq_info, fastq_fwd, fastq_rev]

        }
        .branch { meta_fastq, fastq_info, fastq_fwd, fastq_rev ->
            dna: meta_fastq.sequence_type == 'dna'
            rna: meta_fastq.sequence_type == 'rna'
        }

    //
    // STEP: Handle outputs
    //
    // Set outputs, restoring original meta
    // channel: [ meta, fastq_info, fastq_fwd, fastq_rev ]
    ch_outputs_dna = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_fastq_processed_sorted.dna, ch_inputs),
            ch_inputs_dna_sorted.skip.map { meta -> [meta, [:], [], []] },
        )

    // channel: [ meta, fastq_info, fastq_fwd, fastq_rev ]
    ch_outputs_rna = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_fastq_processed_sorted.rna, ch_inputs),
            ch_inputs_rna_sorted.skip.map { meta -> [meta, [:], [], []] },
        )

    emit:
    fastq_dna = ch_outputs_dna // channel: [ meta, fastq_info, fastq_fwd, fastq_rev ]
    fastq_rna = ch_outputs_rna // channel: [ meta, fastq_info, fastq_fwd, fastq_rev ]
}
