//
// Lift RNA alignments back to genome coordinates
//

include { TARS } from '../../../modules/local/tars/main'

workflow TARS_ALIGNMENT {
    take:
    // Sample data
    ch_inputs           // channel: [mandatory] [ meta ]
    ch_rna_tumor        // channel: [mandatory] [ meta, [aln, ...], [] ]

    // Reference data
    genome_fasta        // channel: [mandatory] /path/to/genome_fasta
    genome_version      // channel: [mandatory] genome version
    genome_fai          // channel: [mandatory] /path/to/genome_fai
    genome_dict         // channel: [mandatory] /path/to/genome_dict
    contigs_mapping_rna // channel: [mandatory] /path/to/contigs_mapping_rna
    unmap_regions_rna   // channel: [mandatory] /path/to/unmap_regions_rna

    main:
    // Sort inputs
    // channel: runnable: [ meta, [aln, ...] ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_rna_tumor
        .branch { meta, alns, _idxs ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.ALN_RNA_TUMOR)
            runnable: alns && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_tars, [aln, ...] ]
    ch_tars_inputs = ch_inputs_sorted.runnable
        .map { meta, alns, _idxs ->

            def sample_id = Utils.getTumorRnaSampleName(meta)

            def meta_tars = [
                key: meta.group_id,
                id: "${meta.group_id}_${sample_id}",
                sample_id: sample_id,
            ]

            return [meta_tars, alns]
        }

    // Run process
    TARS(
        ch_tars_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        contigs_mapping_rna,
        unmap_regions_rna,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, aln, idx ]
    ch_outputs = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('tars_bam'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
        )

    emit:
    rna = ch_outputs // channel: [ meta, aln, idx ]
}
