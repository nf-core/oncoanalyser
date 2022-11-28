//
// LILAC is a WGS tool for HLA typing and somatic CNV and SNV calling
//

include { EXTRACT_AND_INDEX_CONTIG              } from '../../modules/local/custom/lilac_extract_and_index_contig/main'
include { REALIGN_READS                         } from '../../modules/local/custom/lilac_realign_reads_lilac/main'
include { SLICE                                 } from '../../modules/local/custom/lilac_slice/main'

include { LILAC as LILAC_PROCESS                } from '../../modules/local/lilac/main'

workflow LILAC {
  take:
    ch_inputs_bams              // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai]
    ch_inputs_purple_dir        // channel: [val(meta), purple_dir]
    ref_data_genome_fa          //    file: /path/to/genome_fa
    ref_data_genome_fai         //    file: /path/to/genome_fai
    ref_data_genome_version     //     val: genome version
    ref_data_lilac_resource_dir //    file: /path/to/lilac_resource_dir/

  main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Slice HLA regions
    if (ref_data_genome_version == '38') {
      slice_bed = file("${params.ref_data_lilac_resource_dir}/hla.38.alt.umccr.bed", checkIfExists: true)
    } else {
      slice_bed = file("${params.ref_data_lilac_resource_dir}/hla.${ref_data_genome_version}.bed", checkIfExists: true)
    }
    // NOTE(SW): here I remove duplicate files so that we only process each input once
    // NOTE(SW): orphaned reads are sometimes obtained, this is the slicing procedure used
    // in Pipeline5, see LilacBamSlicer.java#L115
    // channel: [val(meta_lilac), bam, bai, bed]
    ch_slice_inputs = WorkflowLilac.get_slice_inputs(ch_inputs_bams, slice_bed)
    // channel: [val(meta_lilac), bam, bai, bed]
    ch_slice_inputs = WorkflowLilac.get_unique_input_files(ch_slice_inputs)
    SLICE(
      ch_slice_inputs,
    )
    ch_versions = ch_versions.mix(SLICE.out.versions)

    // Realign contigs if using hg38 (use of ALT contigs implied)
    if (ref_data_genome_version == '38') {
      // Align reads with chr6
      // NOTE(SW): the aim of this process is to take reads mapping to ALT contigs and align them
      // to the three relevant HLA genes on chr6. All reads including those previously mapped to chr6
      // are realigned for consistency.
      EXTRACT_AND_INDEX_CONTIG(
        'chr6',
        ref_data_genome_fa,
        ref_data_genome_fai,
      )
      REALIGN_READS(
        SLICE.out.bam,
        EXTRACT_AND_INDEX_CONTIG.out.contig,
        EXTRACT_AND_INDEX_CONTIG.out.bwa_indices,
      )

      // Create input channel for LILAC
      // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai, purple_dir]
      ch_lilac_inputs = WorkflowHmftools.group_by_meta(
        WorkflowLilac.sort_slices(REALIGN_READS.out.bam),
        ch_inputs_purple_dir,
      )
    } else {
      // Create input channel for LILAC
      // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai, purple_dir]
      ch_lilac_inputs = WorkflowHmftools.group_by_meta(
        WorkflowLilac.sort_slices(SLICE.out.bam),
        ch_inputs_purple_dir,
      )
    }

    // Run LILAC
    LILAC_PROCESS(
      ch_lilac_inputs,
      ref_data_genome_fa,
      ref_data_genome_version,
      ref_data_lilac_resource_dir,
    )
    ch_versions = ch_versions.mix(LILAC_PROCESS.out.versions)

  emit:
    results = LILAC_PROCESS.out.lilac_dir // channel: [val(meta), lilac_dir]

    versions = ch_versions                // channel: [versions.yml]
}
