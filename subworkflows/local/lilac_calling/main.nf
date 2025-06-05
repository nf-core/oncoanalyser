//
// LILAC is a WGS tool for HLA typing and somatic CNV and SNV calling
//

import Constants
import Utils

include { CUSTOM_EXTRACTCONTIG as EXTRACTCONTIG } from '../../../modules/local/custom/lilac_extract_and_index_contig/main'
include { CUSTOM_REALIGNREADS as REALIGNREADS   } from '../../../modules/local/custom/lilac_realign_reads_lilac/main'
include { CUSTOM_SLICE as SLICEBAM              } from '../../../modules/local/custom/lilac_slice/main'
include { LILAC                                 } from '../../../modules/local/lilac/main'

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

    // Realign reads mapping to HLA regions and homologus regions if using reference genome with ALT contigs
    // NOTE(SW): the aim of this process is to take reads mapping to ALT contigs and align them to the three
    // relevant HLA genes on chr6. All reads including those previously mapped to chr6 are realigned for
    // consistency.
    if (params.genome_type == 'alt') {

        // Flatten into BAM/BAI pairs, select inputs that are eligible to run
        // channel: runnable: [ meta_extra, bam, bai ]
        // channel: skip: [ meta_extra ]
        ch_realign_inputs_sorted = ch_dna_inputs_sorted.runnable
            .flatMap { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->

                def tumor_sample_id = Utils.hasTumorDna(meta) ? Utils.getTumorDnaSampleName(meta) : []
                def normal_sample_id = Utils.hasNormalDna(meta) ? Utils.getNormalDnaSampleName(meta) : []

                return [
                    [[key: meta.group_id, *:meta, sample_id: tumor_sample_id, sample_type: 'tumor'], tumor_bam, tumor_bai],
                    [[key: meta.group_id, *:meta, sample_id: normal_sample_id, sample_type: 'normal'], normal_bam, normal_bai],
                ]
            }
            .branch { meta_extra, bam, bai ->
                runnable: bam && bai
                skip: true
                    return meta_extra
            }

        //
        // MODULE: Custom BAM slice (LILAC)
        //
        // Create process input channel
        // channel: [ meta_realign, bam, bai ]
        ch_slice_inputs = ch_realign_inputs_sorted.runnable
            .map { meta_extra, bam, bai ->

                def meta_realign = [
                    key: meta_extra.group_id,
                    id: "${meta_extra.group_id}__${meta_extra.sample_id}",
                    sample_id: meta_extra.sample_id,
                    sample_type: meta_extra.sample_type,
                ]

                return [meta_realign, bam, bai]
            }

        // Run process
        SLICEBAM(
            ch_slice_inputs,
            hla_slice_bed,
        )

        ch_versions = ch_versions.mix(SLICEBAM.out.versions)

        //
        // MODULE: Custom extract contig (LILAC)
        //
        // Only run if we have runnable inputs, no blocking since operating only on input metas
        ch_extract_contig_run = ch_realign_inputs_sorted.runnable
            .toList()
            .map { !it.isEmpty() }

        EXTRACTCONTIG(
            'chr6',
            genome_fasta,
            genome_fai,
            ch_extract_contig_run,
        )

        ch_versions = ch_versions.mix(EXTRACTCONTIG.out.versions)

        //
        // MODULE: Custom realign reads (LILAC)
        //
        REALIGNREADS(
            SLICEBAM.out.bam,
            EXTRACTCONTIG.out.contig,
            EXTRACTCONTIG.out.bwamem2_index,
        )

        ch_versions = ch_versions.mix(REALIGNREADS.out.versions)

        // Separate all BAMs by sample type so they can be merged with desired order
        // channel: [ < meta_extra OR meta_realign >, bam, bai ]
        ch_slice_reunited_bams = Channel.empty()
            .mix(
                ch_realign_inputs_sorted.skip.map { meta_extra -> [meta_extra, [], []] },
                REALIGNREADS.out.bam,
            )
            .branch { meta_ambiguous, bam, bai ->
                tumor: meta_ambiguous.sample_type == 'tumor'
                normal: meta_ambiguous.sample_type == 'normal'
            }

        // Restore meta, pair tumor and normal BAMs
        // channel: [ meta, tumor_dna_bam, tumor_dna_bai, normal_dna_bam, normal_dna_bai ]
        ch_dna_inputs_ready = WorkflowOncoanalyser.groupByMeta(
            WorkflowOncoanalyser.restoreMeta(ch_slice_reunited_bams.tumor, ch_inputs),
            WorkflowOncoanalyser.restoreMeta(ch_slice_reunited_bams.normal, ch_inputs),
        )

    } else {

        // channel: [ meta, tumor_dna_bam, tumor_dna_bai, normal_dna_bam, normal_dna_bai ]
        ch_dna_inputs_ready = ch_dna_inputs_sorted.runnable

    }

    //
    // MODULE: LILAC
    //
    // Create process input channel
    // channel: [ meta_lilac, normal_dna_bam, normal_dna_bai, tumor_dna_bam, tumor_dna_bai, tumor_rna_bam, tumor_rna_bai, purple_dir ]
    ch_lilac_inputs = WorkflowOncoanalyser.groupByMeta(
        ch_dna_inputs_ready,
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
