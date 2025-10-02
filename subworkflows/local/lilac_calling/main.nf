//
// LILAC is a WGS tool for HLA typing and somatic CNV and SNV calling
//

import Constants
import Utils

include { LILAC } from '../../../modules/local/lilac/main'

workflow LILAC_CALLING {
    take:
    // Sample data
    ch_inputs          // channel: [mandatory] [ meta ]
    ch_tumor_bam       // channel: [mandatory] [ meta, bam, bai ]
    ch_normal_bam      // channel: [mandatory] [ meta, bam, bai ]
    ch_tumor_rna_bam   // channel: [mandatory] [ meta, bam, bai ]
    ch_purple          // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_fasta       // channel: [mandatory] /path/to/genome_fasta
    genome_version     // channel: [mandatory] genome version
    genome_fai         // channel: [mandatory] /path/to/genome_fai
    lilac_resource_dir // channel: [mandatory] /path/to/lilac_resource_dir/
    hla_slice_bed      // channel: [mandatory] /path/to/hla_slice_bed
    targeted_mode   // boolean: [mandatory] Running in targeted/panel mode?

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources and sort for DNA BAMs
    // channel: runnable: [ meta, tumor_dna_bam, tumor_dna_bai, normal_dna_bam, normal_dna_bai ]
    // channel: skip: [ meta ]
    ch_dna_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_tumor_bam,
        ch_normal_bam,
    )
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            return [
                meta,
                Utils.selectCurrentOrExisting(tumor_bam, meta, Constants.INPUT.BAM_REDUX_DNA_TUMOR),
                tumor_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_TUMOR),
                Utils.selectCurrentOrExisting(normal_bam, meta, Constants.INPUT.BAM_REDUX_DNA_NORMAL),
                normal_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_NORMAL),
            ]
        }
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.LILAC_DIR)

            runnable: (tumor_bam || normal_bam) && !has_existing
            skip: true
                return meta
        }

    //
    // MODULE: LILAC
    //
    // Create process input channel
    // channel: [ meta_lilac, normal_dna_bam, normal_dna_bai, tumor_dna_bam, tumor_dna_bai, tumor_rna_bam, tumor_rna_bai, purple_dir ]
    ch_lilac_inputs = WorkflowOncoanalyser.groupByMeta(
        ch_dna_inputs_sorted.runnable,
        ch_tumor_rna_bam,
        ch_purple,
    )
        .map { meta, tbam_dna, tbai_dna, nbam_dna, nbai_dna, tbam_rna, tbai_rna, purple_dir ->

            def meta_lilac = [
                key: meta.group_id,
                id: meta.group_id,
            ]

            if (Utils.hasTumorDna(meta)) {
                meta_lilac.tumor_id = Utils.getTumorDnaSampleName(meta)
            }

            if (Utils.hasNormalDna(meta)) {
                meta_lilac.normal_id = Utils.getNormalDnaSampleName(meta)
            }

            return [
                meta_lilac,
                nbam_dna,
                nbai_dna,
                tbam_dna,
                tbai_dna,
                Utils.selectCurrentOrExisting(tbam_rna, meta, Constants.INPUT.BAM_RNA_TUMOR),
                Utils.selectCurrentOrExisting(tbai_rna, meta, Constants.INPUT.BAI_RNA_TUMOR),
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
            ]
        }

    // Run process
    LILAC(
        ch_lilac_inputs,
        genome_fasta,
        genome_fai,
        genome_version,
        lilac_resource_dir,
        targeted_mode
    )

    ch_versions = ch_versions.mix(LILAC.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, amber_dir ]
    ch_outputs = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(LILAC.out.lilac_dir, ch_inputs),
            ch_dna_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    lilac_dir = ch_outputs  // channel: [ meta, lilac_dir ]

    versions  = ch_versions // channel: [ versions.yml ]
}
