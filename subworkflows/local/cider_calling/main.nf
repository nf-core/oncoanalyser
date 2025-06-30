//
// CIDER identifies and annotates CDR3 sequences of IG and TCR loci
//

import Constants
import Utils

include { CIDER } from '../../../modules/local/cider/main'

workflow CIDER_CALLING {
    take:
    // Sample data
    ch_inputs        // channel: [mandatory] [ meta ]
    ch_tumor_dna_bam // channel: [mandatory] [ meta, bam, bai ]
    ch_tumor_rna_bam // channel: [mandatory] [ meta, bam, bai ]

    // Reference data
    genome_version // channel: [mandatory] genome version
    human_blastdb  // channel: [mandatory] /path/to/human_blastdb

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Sort inputs, separate by DNA and RNA
    // channel: runnable: [ meta, bam, bai ]
    // channel: skip: [ meta ]
    ch_inputs_tumor_dna_sorted = ch_tumor_dna_bam
        .map { meta, bam, bai ->
            return [
                meta,
                Utils.selectCurrentOrExisting(bam, meta, Constants.INPUT.BAM_REDUX_DNA_TUMOR),
                bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_TUMOR),
            ]
        }
        .branch { meta, bam, bai ->
            runnable: bam
            skip: true
                return meta
        }

    // channel: runnable: [ meta, bam, bai ]
    // channel: skip: [ meta ]
    ch_inputs_tumor_rna_sorted = ch_tumor_rna_bam
        .map { meta, bam, bai ->
            return [
                meta,
                Utils.selectCurrentOrExisting(bam, meta, Constants.INPUT.BAM_RNA_TUMOR),
                bai ?: Utils.getInput(meta, Constants.INPUT.BAI_RNA_TUMOR),
            ]
        }
        .branch { meta, bam, bai ->
            runnable: bam
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_cider, bam, bai ]
    ch_cider_inputs = Channel.empty()
        .mix(
            ch_inputs_tumor_dna_sorted.runnable.map { meta, bam, bai -> [meta, Utils.getTumorDnaSample(meta), bam, bai] },
            ch_inputs_tumor_rna_sorted.runnable.map { meta, bam, bai -> [meta, Utils.getTumorRnaSample(meta), bam, bai] },
        )
        .map { meta, meta_sample, bam, bai ->

            def meta_cider = [
                key: meta.group_id,
                id: "${meta.group_id}_${meta_sample.sample_id}",
                sample_id: meta_sample.sample_id,
            ]

            return [meta_cider, bam, bai]
        }

    // Run process
    CIDER(
        ch_cider_inputs,
        genome_version,
        human_blastdb,
    )

    ch_versions = ch_versions.mix(CIDER.out.versions)

    emit:
    versions  = ch_versions // channel: [ versions.yml ]
}
