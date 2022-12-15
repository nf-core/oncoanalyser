//
// LILAC is a WGS tool for HLA typing and somatic CNV and SNV calling
//
import Constants

include { EXTRACT_AND_INDEX_CONTIG              } from '../../modules/local/custom/lilac_extract_and_index_contig/main'
include { REALIGN_READS                         } from '../../modules/local/custom/lilac_realign_reads_lilac/main'
include { SLICE                                 } from '../../modules/local/custom/lilac_slice/main'

include { LILAC as LILAC_PROCESS                } from '../../modules/local/lilac/main'

workflow LILAC {
    take:
        ch_inputs_bams              // channel: [val(meta_lilac), tumor_bam, normal_bam, tumor_bai, normal_bai]
        ch_inputs_purple_dir        // channel: [val(meta_lilac), purple_dir]
        ref_data_genome_fasta       //    file: /path/to/genome_fasta
        ref_data_genome_fai         //    file: /path/to/genome_fai
        ref_data_lilac_resource_dir //    file: /path/to/lilac_resource_dir/

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Slice HLA regions
        ch_slice_bed = ref_data_lilac_resource_dir
            .map { lilac_dir ->
                def filename
                if (params.ref_data_genome_type == 'no_alt') {
                    filename = "hla.${params.ref_data_genome_version}.bed"
                } else {
                    filename = "hla.38.alt.umccr.bed"
                }

                def filepath = lilac_dir.resolve(filename).toUri().toString()
                // NOTE(SW): S3Path.toUri() includes an extra forward slash in the URI or path,
                // which is removed here. Unsure whether this occurs in other URIs types
                if (filepath.startsWith('s3:///')) {
                    filepath = filepath.replaceFirst(/^s3:\/\/\//, 's3://')
                }
                return file(filepath, checkIfExists: true)
            }
        // NOTE(SW): here I remove duplicate files so that we only process each input once
        // NOTE(SW): orphaned reads are sometimes obtained, this is the slicing procedure used
        // in Pipeline5, see LilacBamSlicer.java#L115
        // channel: [val(meta_lilac), bam, bai, bed]
        ch_slice_inputs = WorkflowLilac.getSliceInputs(ch_inputs_bams, ch_slice_bed)
        // channel: [val(meta_lilac), bam, bai, bed]
        ch_slice_inputs = WorkflowLilac.getUniqueInputFiles(ch_slice_inputs)
        SLICE(
            ch_slice_inputs,
        )
        ch_versions = ch_versions.mix(SLICE.out.versions)

        // Realign contigs if using 38 (use of ALT contigs implied)
        // channel: [val(meta)]
        ch_metas = ch_inputs_bams.map { return it[0] }
        if (params.ref_data_genome_type == 'alt') {
            // Align reads with chr6
            // NOTE(SW): the aim of this process is to take reads mapping to ALT contigs and align them
            // to the three relevant HLA genes on chr6. All reads including those previously mapped to chr6
            // are realigned for consistency.
            EXTRACT_AND_INDEX_CONTIG(
                'chr6',
                ref_data_genome_fasta,
                ref_data_genome_fai,
            )
            REALIGN_READS(
                SLICE.out.bam,
                EXTRACT_AND_INDEX_CONTIG.out.contig,
                EXTRACT_AND_INDEX_CONTIG.out.bwa_indices,
            )

            // Create input channel for LILAC
            // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai]
            ch_lilac_inputs_slices = WorkflowOncoanalyser.restoreMeta(
                WorkflowLilac.sortSlices(REALIGN_READS.out.bam),
                ch_metas,
            )
            // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai, purple_dir]
            ch_lilac_inputs_full = WorkflowOncoanalyser.groupByMeta(
                ch_lilac_inputs_slices,
                ch_inputs_purple_dir,
            )
        } else {
            // Create input channel for LILAC
            // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai]
            ch_lilac_inputs_slices = WorkflowOncoanalyser.restoreMeta(
                WorkflowLilac.sortSlices(SLICE.out.bam),
                ch_metas,
            )
            // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai, purple_dir]
            ch_lilac_inputs_full = WorkflowOncoanalyser.groupByMeta(
                ch_lilac_inputs_slices,
                ch_inputs_purple_dir,
            )
        }

        // Run LILAC
        // channel: [val(meta_lilac), tumor_bam, normal_bam, tumor_bai, normal_bai, purple_dir]
        ch_lilac_inputs = ch_lilac_inputs_full
            .map {
                def meta = it[0]
                def meta_lilac = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: meta.get(['sample_name', Constants.DataType.TUMOR]),
                    normal_id: meta.get(['sample_name', Constants.DataType.NORMAL]),
                ]
                return [meta_lilac, *it[1..-1]]
            }
        LILAC_PROCESS(
            ch_lilac_inputs,
            ref_data_genome_fasta,
            params.ref_data_genome_version,
            ref_data_lilac_resource_dir,
        )
        ch_versions = ch_versions.mix(LILAC_PROCESS.out.versions)

    emit:
        results = LILAC_PROCESS.out.lilac_dir // channel: [val(meta), lilac_dir]

        versions = ch_versions                // channel: [versions.yml]
}
