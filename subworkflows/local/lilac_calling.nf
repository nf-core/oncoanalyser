//
// LILAC is a WGS tool for HLA typing and somatic CNV and SNV calling
//
import Constants
import Utils

include { CUSTOM_EXTRACTCONTIG as EXTRACTCONTIG } from '../../modules/local/custom/lilac_extract_and_index_contig/main'
include { CUSTOM_REALIGNREADS as REALIGNREADS   } from '../../modules/local/custom/lilac_realign_reads_lilac/main'
include { CUSTOM_SLICE as SLICEBAM              } from '../../modules/local/custom/lilac_slice/main'
include { LILAC                                 } from '../../modules/local/lilac/main'

include { CHANNEL_GROUP_INPUTS } from './channel_group_inputs'

workflow LILAC_CALLING {
    take:
        // Sample data
        ch_inputs
        ch_purple

        // Reference data
        ref_data_genome_fasta       //    file: /path/to/genome_fasta
        ref_data_genome_fai         //    file: /path/to/genome_fai
        ref_data_lilac_resource_dir //    file: /path/to/lilac_resource_dir/
        ref_data_hla_slice_bed

        // Params
        run

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Get input meta groups
        CHANNEL_GROUP_INPUTS(
            ch_inputs,
        )

        // Select input sources
        // channel: [val(meta), purple_dir]
        if (run.purple) {
            ch_lilac_inputs_purple = Channel.empty()
                .mix(
                    ch_purple,
                    CHANNEL_GROUP_INPUTS.out.wgs_absent.map { meta -> [meta, []] },
                )
        } else {
            ch_lilac_inputs_purple = WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PURPLE_DIR, type: 'optional')
        }

        // Create channel with all available input BAMs
        // First obtain WTS BAMs
        // channel: [val(meta), wts_bam, wts_bai]
        ch_lilac_bams_wts = Channel.empty()
            .mix(
                CHANNEL_GROUP_INPUTS.out.wts_present.map { meta ->
                    def bam = Utils.getTumorWtsBam(meta)
                    return [meta, bam, "${bam}.bai"]
                },
                CHANNEL_GROUP_INPUTS.out.wts_absent.map { meta -> [meta, [], []] },
            )

        ch_lilac_bams_wgs = Channel.empty()
            .mix(
                CHANNEL_GROUP_INPUTS.out.wgs_present.map { meta ->
                    def tumor_bam = Utils.getTumorWgsBam(meta)
                    def normal_bam = Utils.getNormalWgsBam(meta)
                    [meta, tumor_bam, normal_bam, "${tumor_bam}.bai", "${normal_bam}.bai"]
                },
                CHANNEL_GROUP_INPUTS.out.wgs_absent.map { meta -> [meta, [], [], [], []] },
            )

        // Combine WGS and WTS BAMs
        // channel: [val(meta), normal_wgs_bam, normal_wgs_bai, tumor_wgs_bam, tumor_wgs_bai, tumor_wts_bam, tumor_wts_bai]
        ch_lilac_bams = WorkflowOncoanalyser.groupByMeta(
            ch_lilac_bams_wts,
            ch_lilac_bams_wgs,
            flatten_mode: 'nonrecursive',
        )
            .map { data ->
                def meta = data[0]
                def (tbam_wts, tbai_wts, tbam_wgs, nbam_wgs, tbai_wgs, nbai_wgs) = data[1..-1]
                return [meta, nbam_wgs, nbai_wgs, tbam_wgs, tbai_wgs, tbam_wts, tbai_wts]
            }

        // Slice HLA regions
        // NOTE(SW): here I remove duplicate files so that we only process each input once
        // NOTE(SW): orphaned reads are sometimes obtained, this is the slicing procedure used
        // in Pipeline5, see LilacBamSlicer.java#L115
        // channel: [val(meta_lilac), bam, bai]
        ch_slice_inputs = WorkflowLilac.getSliceInputs(ch_lilac_bams)
        // Isolate meta containing expected file count to use for non-blocking groupTuple later
        ch_slice_meta_individual = ch_slice_inputs
            .map {
                def meta_lilac = it[0]
                return [key: meta_lilac.key, count: meta_lilac.count]
            }
        // Apply slicing to unique files only
        // channel: [val(meta_lilac), bam, bai]
        ch_slice_inputs_unique = WorkflowLilac.getUniqueInputFiles(ch_slice_inputs)
        SLICEBAM(
            ch_slice_inputs_unique,
            ref_data_hla_slice_bed,
        )
        ch_versions = ch_versions.mix(SLICEBAM.out.versions)

        // Realign contigs if using 38 (use of ALT contigs implied)
        // channel: [val(meta_lilac), bam, bai]
        ch_slices_out = SLICEBAM.out.bam
        if (params.ref_data_genome_type == 'alt') {
            // Align reads with chr6
            // NOTE(SW): the aim of this process is to take reads mapping to ALT contigs and align them
            // to the three relevant HLA genes on chr6. All reads including those previously mapped to chr6
            // are realigned for consistency.
            EXTRACTCONTIG(
                'chr6',
                ref_data_genome_fasta,
                ref_data_genome_fai,
            )
            ch_versions = ch_versions.mix(EXTRACTCONTIG.out.versions)

            ch_slice_bams = SLICEBAM.out.bam
                .branch { meta_lilac, bam, bai ->
                    def sequence_type = Utils.getEnumFromString(meta_lilac.sequence_type_str, Constants.SequenceType)
                    wgs: sequence_type == Constants.SequenceType.WGS
                    wts: sequence_type == Constants.SequenceType.WTS
                }

            REALIGNREADS(
                ch_slice_bams.wgs,
                EXTRACTCONTIG.out.contig,
                EXTRACTCONTIG.out.bwa_indices,
            )
            ch_slices_out = Channel.empty().mix(REALIGNREADS.out.bam, ch_slice_bams.wts)
            ch_versions = ch_versions.mix(REALIGNREADS.out.versions)
        }

        // Re-replicate and flow expected file count into meta
        // channel: [val(meta_lilac), [sequence_type_str, sample_type_str, bam, bai]]
        ch_slices_out_individual = ch_slices_out
            .flatMap{ data ->
                def meta_lilac = data[0]
                def fps = data[1..-1]
                meta_lilac.keys.collect { key ->
                    return [[key: key], [meta_lilac.sequence_type_str, meta_lilac.sample_type_str, *fps]]
                }
            }
        // Adding expected file count
        // channel: [val(meta_lilac), [sequence_type_str, sample_type_str, bam, bai]]
        ch_slices_ready = WorkflowOncoanalyser.joinMeta(
            ch_slices_out_individual,
            ch_slice_meta_individual,
            key_a: 'key',
        )

        // Gather and order files from same grouping using non-blocked groupTuple via provided file counts
        // channel: [val(meta_lilac), normal_wgs_bam, normal_wgs_bai, tumor_wgs_bam, tumor_wgs_bai, tumor_wts_bam, tumor_wts_bai]
        ch_slices_organised = WorkflowLilac.sortSlices(ch_slices_ready)

        // Restore original meta so we can join with PURPLE directory
        // channel: [val(meta)]
        ch_metas = ch_lilac_bams.map { return it[0] }
        // channel: [val(meta), normal_wgs_bam, normal_wgs_bai, tumor_wgs_bam, tumor_wgs_bai, tumor_wts_bam, tumor_wts_bai]
        ch_lilac_inputs_slices = WorkflowOncoanalyser.restoreMeta(
            ch_slices_organised,
            ch_metas,
        )

        // Get inputs from PURPLE
        // channel: [val(meta), gene_cn]
        ch_lilac_inputs_gene_cn = ch_purple
            .map { meta, purple_dir ->
                if (purple_dir == []) {
                    return [meta, []]
                }

                def tumor_id = Utils.getTumorWgsSampleName(meta)
                def gene_cn = file(purple_dir).resolve("${tumor_id}.purple.cnv.gene.tsv")
                return gene_cn.exists() ? [meta, gene_cn] : [meta, []]
            }

        // channel: [val(meta), smlv_vcf]
        ch_lilac_inputs_smlv_vcf = ch_purple
            .map { meta, purple_dir ->
                if (purple_dir == []) {
                    return [meta, []]
                }

                def tumor_id = Utils.getTumorWgsSampleName(meta)
                def smlv_vcf = file(purple_dir).resolve("${tumor_id}.purple.somatic.vcf.gz")
                return smlv_vcf.exists() ? [meta, smlv_vcf] : [meta, []]
            }

        // channel: [val(meta), normal_wgs_bam, normal_wgs_bai, tumor_wgs_bam, tumor_wgs_bai, tumor_wts_bam, tumor_wts_bai, purple_dir]
        ch_lilac_inputs_full = WorkflowOncoanalyser.groupByMeta(
            ch_lilac_inputs_slices,
            ch_lilac_inputs_gene_cn,
            ch_lilac_inputs_smlv_vcf,
            flatten_mode: 'nonrecursive',
        )

        // Create final input channel for LILAC, remove samples with only WTS BAMs
        // channel: [val(meta_lilac), normal_wgs_bam, normal_wgs_bai, tumor_wgs_bam, tumor_wgs_bai, tumor_wts_bam, tumor_wts_bai, gene_cn, smlv_vcf]]
        ch_lilac_inputs = ch_lilac_inputs_full
            .map {
                def meta = it[0]
                def fps = it[1..-1]

                // LILAC requires either tumor or normal WGS BAM
                if (fps[0] == [] && fps[2] == []) {
                    return Constants.META_PLACEHOLDER
                }

                def tumor_id_key = ['sample_name', Constants.SampleType.TUMOR, Constants.SequenceType.WGS]
                def normal_id_key = ['sample_name', Constants.SampleType.NORMAL, Constants.SequenceType.WGS]

                def meta_lilac = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: meta.containsKey(tumor_id_key) ? meta.getAt(tumor_id_key) : '',
                    normal_id: meta.containsKey(normal_id_key) ? meta.getAt(normal_id_key) : '',
                ]
                return [meta_lilac, *fps]
            }
            .filter { it != Constants.META_PLACEHOLDER }

        // Run LILAC
        LILAC(
            ch_lilac_inputs,
            ref_data_genome_fasta,
            params.ref_data_genome_version,
            ref_data_lilac_resource_dir,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(LILAC.out.versions)
        ch_outputs = WorkflowOncoanalyser.restoreMeta(LILAC.out.lilac_dir, ch_inputs)

    emit:
        lilac_dir = ch_outputs // channel: [val(meta), lilac_dir]

        versions = ch_versions // channel: [versions.yml]
}
