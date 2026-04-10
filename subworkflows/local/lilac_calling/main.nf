//
// LILAC is a WGS tool for HLA typing and somatic CNV and SNV calling
//

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
    ch_dna_inputs_sorted = channels.WorkflowChannels.groupByMeta(
        ch_tumor_bam,
        ch_normal_bam,
    )
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            return [
                meta,
                sample.Inputs.preferUserProvidedInput(tumor_bam, meta, sample.FileKey.BAM_REDUX_DNA_TUMOR),
                sample.Inputs.preferPipelineOutput(tumor_bai, meta, sample.FileKey.BAI_DNA_TUMOR),
                sample.Inputs.preferUserProvidedInput(normal_bam, meta, sample.FileKey.BAM_REDUX_DNA_NORMAL),
                sample.Inputs.preferPipelineOutput(normal_bai, meta, sample.FileKey.BAI_DNA_NORMAL),
            ]
        }
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->

            def has_existing = sample.Inputs.hasExisting(meta, sample.FileKey.LILAC_DIR)

            runnable: (tumor_bam || normal_bam) && !has_existing
            skip: true
                return meta
        }

    //
    // MODULE: LILAC
    //
    // Create process input channel
    // channel: [ meta_lilac, normal_dna_bam, normal_dna_bai, tumor_dna_bam, tumor_dna_bai, tumor_rna_bam, tumor_rna_bai, purple_dir ]
    ch_lilac_inputs = channels.WorkflowChannels.groupByMeta(
        ch_dna_inputs_sorted.runnable,
        ch_tumor_rna_bam,
        ch_purple,
    )
        .map { meta, tbam_dna, tbai_dna, nbam_dna, nbai_dna, tbam_rna, tbai_rna, purple_dir ->

            def meta_lilac = [
                key: meta.group_id,
                id: meta.group_id,
            ]

            if (sample.Inputs.hasTumorDna(meta)) {
                meta_lilac.tumor_id = sample.Inputs.getTumorDnaSampleName(meta)
            }

            if (sample.Inputs.hasNormalDna(meta)) {
                meta_lilac.normal_id = sample.Inputs.getNormalDnaSampleName(meta)
            }

            return [
                meta_lilac,
                nbam_dna,
                nbai_dna,
                tbam_dna,
                tbai_dna,
                sample.Inputs.preferUserProvidedInput(tbam_rna, meta, sample.FileKey.BAM_RNA_TUMOR),
                sample.Inputs.preferUserProvidedInput(tbai_rna, meta, sample.FileKey.BAI_RNA_TUMOR),
                sample.Inputs.preferUserProvidedInput(purple_dir, meta, sample.FileKey.PURPLE_DIR),
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
            channels.WorkflowChannels.restoreMeta(LILAC.out.lilac_dir, ch_inputs),
            channels.PlaceholderChannels.toolDir(ch_dna_inputs_sorted.skip),
        )

    emit:
    lilac_dir = ch_outputs  // channel: [ meta, lilac_dir ]

    versions  = ch_versions // channel: [ versions.yml ]
}
