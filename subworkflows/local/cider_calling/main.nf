//
// CIDER identifies and annotates CDR3 sequences of IG and TCR loci
//

include { CIDER } from '../../../modules/local/cider/main'

workflow CIDER_CALLING {
    take:
    // Sample data
    ch_inputs          // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor // channel: [mandatory] [ meta, redux_dir ]
    ch_tumor_rna_aln   // channel: [mandatory] [ meta, aln, idx ]

    // Reference data
    genome_fasta       // channel: [mandatory] /path/to/genome_fasta
    genome_version     // channel: [mandatory] genome version
    genome_fai         // channel: [mandatory] /path/to/genome_fai
    genome_dict        // channel: [mandatory] /path/to/genome_dict
    genome_img         // channel: [optional]  /path/to/genome_img

    main:
    // Select input sources then sort, separate by DNA and RNA
    // channel: runnable: [ meta, aln, idx ]
    // channel: skip: [ meta ]
    ch_inputs_tumor_dna_sorted = ch_redux_dir_tumor
        .map { meta, redux_dir_tumor ->

            def redux_dir_tumor_selected = Utils.selectCurrentOrExisting(redux_dir_tumor, meta, Constants.INPUT.REDUX_DIR_TUMOR)
            def (tumor_aln, tumor_idx) = Utils.getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)

            return [meta, tumor_aln, tumor_idx]

        }
        .branch { meta, aln, idx ->
            runnable: aln
            skip: true
                return meta
        }

    // channel: runnable: [ meta, aln, idx ]
    // channel: skip: [ meta ]
    ch_inputs_tumor_rna_sorted = ch_tumor_rna_aln
        .map { meta, aln, idx ->
            return [
                meta,
                Utils.selectCurrentOrExisting(aln, meta, Constants.INPUT.ALN_RNA_TUMOR),
                idx ?: Utils.getInput(meta, Constants.INPUT.IDX_RNA_TUMOR),
            ]
        }
        .branch { meta, aln, idx ->
            runnable: aln
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_cider, aln, idx ]
    ch_cider_inputs = channel.empty()
        .mix(
            ch_inputs_tumor_dna_sorted.runnable.map { meta, aln, idx -> [meta, Utils.getTumorDnaSample(meta), aln, idx] },
            ch_inputs_tumor_rna_sorted.runnable.map { meta, aln, idx -> [meta, Utils.getTumorRnaSample(meta), aln, idx] },
        )
        .map { meta, meta_sample, aln, idx->

            def meta_cider = [
                key: meta.group_id,
                id: "${meta.group_id}_${meta_sample.sample_id}",
                sample_id: meta_sample.sample_id,
            ]

            return [meta_cider, aln, idx]
        }

    // Run process
    CIDER(
        ch_cider_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        genome_img,
    )
}
