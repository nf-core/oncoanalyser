//
// LILAC is a WGS tool for HLA typing and somatic CNV and SNV calling
//

import Constants
import Inputs

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
    targeted_mode      // boolean: [mandatory] Set targeted mode

    // Params
    sequencing_type    // string:  [mandatory] sequencing type

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
                Inputs.overrideWithExistingInput(tumor_bam, meta, Constants.INPUT.BAM_REDUX_DNA_TUMOR),
                Inputs.fallbackToExistingInput(tumor_bai, meta, Constants.INPUT.BAI_DNA_TUMOR),
                Inputs.overrideWithExistingInput(normal_bam, meta, Constants.INPUT.BAM_REDUX_DNA_NORMAL),
                Inputs.fallbackToExistingInput(normal_bai, meta, Constants.INPUT.BAI_DNA_NORMAL),
            ]
        }
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->

            def has_existing = Inputs.hasExistingInput(meta, Constants.INPUT.LILAC_DIR)

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

            if (Inputs.hasTumorDna(meta)) {
                meta_lilac.tumor_id = Inputs.getTumorDnaSampleName(meta)
            }

            if (Inputs.hasNormalDna(meta)) {
                meta_lilac.normal_id = Inputs.getNormalDnaSampleName(meta)
            }

            return [
                meta_lilac,
                nbam_dna,
                nbai_dna,
                tbam_dna,
                tbai_dna,
                Inputs.overrideWithExistingInput(tbam_rna, meta, Constants.INPUT.BAM_RNA_TUMOR),
                Inputs.overrideWithExistingInput(tbai_rna, meta, Constants.INPUT.BAI_RNA_TUMOR),
                Inputs.overrideWithExistingInput(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
            ]
        }

    // Run process
    LILAC(
        ch_lilac_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        lilac_resource_dir,
        targeted_mode,
        sequencing_type,
    )

    ch_versions = ch_versions.mix(LILAC.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, amber_dir ]
    ch_outputs = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(LILAC.out.lilac_dir, ch_inputs),
            PlaceholderChannels.toolDir(ch_dna_inputs_sorted.skip),
        )

    emit:
    lilac_dir = ch_outputs  // channel: [ meta, lilac_dir ]

    versions  = ch_versions // channel: [ versions.yml ]
}
