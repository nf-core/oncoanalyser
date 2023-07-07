//
// LILAC is a WGS tool for HLA typing and somatic CNV and SNV calling
//

import Constants
import Utils

include { CUSTOM_EXTRACTCONTIG as EXTRACTCONTIG } from '../../modules/local/custom/lilac_extract_and_index_contig/main'
include { CUSTOM_REALIGNREADS as REALIGNREADS   } from '../../modules/local/custom/lilac_realign_reads_lilac/main'
include { CUSTOM_SLICE as SLICEBAM              } from '../../modules/local/custom/lilac_slice/main'
include { LILAC                                 } from '../../modules/local/lilac/main'

workflow LILAC_CALLING {
    take:
        // Sample data
        ch_inputs          // channel: [mandatory] [ meta ]
        ch_purple          // channel: [optional]  [ meta, purple_dir ]

        // Reference data
        genome_fasta       // channel: [mandatory] /path/to/genome_fasta
        genome_version     // channel: [mandatory] genome version
        genome_fai         // channel: [mandatory] /path/to/genome_fai
        lilac_resource_dir // channel: [mandatory] /path/to/lilac_resource_dir/
        hla_slice_bed      // channel: [mandatory] /path/to/hla_slice_bed

        // Params
        run_config         // channel: [mandatory] run configuration

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [ meta, purple_dir ]
        if (run_config.stages.purple) {
            ch_lilac_inputs_purple = ch_purple
        } else {
            ch_lilac_inputs_purple = WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PURPLE_DIR, type: 'optional')
        }

        // Create channel with all available input BAMs
        // channel: [ meta, wts_bam, wts_bai ]
        ch_lilac_bams_wts = Channel.empty()
        if (run_config.mode == Constants.RunMode.WTS || run_config.mode == Constants.RunMode.WGTS) {
            ch_lilac_bams_wts = ch_inputs
                .map { meta ->
                    def bam = Utils.getTumorWtsBam(meta)
                    return [meta, bam, "${bam}.bai"]
                }
        } else {
            ch_lilac_bams_wts = ch_inputs.map { meta -> [meta, [], []] }
        }

        // channel: [ meta, bam, bai ]
        ch_lilac_bams = Channel.empty()
        if (run_config.mode == Constants.RunMode.WGS || run_config.mode == Constants.RunMode.WGTS) {

            ch_lilac_bams = ch_inputs
                .map { meta ->

                    def tumor_bam = Utils.getTumorWgsBam(meta)

                    def normal_bam = []
                    def normal_bai = []

                    if (run_config.type == Constants.RunType.TUMOR_NORMAL) {
                        normal_bam = Utils.getNormalWgsBam(meta)
                        normal_bai = "${normal_bam}.bai"
                    }

                    [meta, tumor_bam, normal_bam, "${tumor_bam}.bai", normal_bai]
                }

        } else {
            ch_lilac_bams = ch_inputs.map { meta -> [meta, [], [], [], []] }
        }

        // Combine BAMs
        // channel: [ meta, normal_wgs_bam, normal_wgs_bai, tumor_bam, tumor_bai, tumor_wts_bam, tumor_wts_bai ]
        ch_lilac_bams = WorkflowOncoanalyser.groupByMeta(
            ch_lilac_bams_wts,
            ch_lilac_bams,
            flatten_mode: 'nonrecursive',
        )
            .map { data ->
                def meta = data[0]
                def (tbam_wts, tbai_wts, tbam, nbam, tbai, nbai) = data[1..-1]
                return [meta, nbam, nbai, tbam, tbai, tbam_wts, tbai_wts]
            }

        // Set tumor sequence type for non-WTS input (i.e. WGS or panel)
        switch(run_config.mode) {
            case Constants.RunMode.WGS:
            case Constants.RunMode.WGTS:
                tumor_sequence_type = Constants.SequenceType.WGS
                break
            default:
                assert false
        }

        // Slice HLA regions
        // NOTE(SW): here I remove duplicate files so that we only process each input once
        // NOTE(SW): orphaned reads are sometimes obtained, this is the slicing procedure used
        // in Pipeline5, see LilacBamSlicer.java#L115
        // channel: [ meta_lilac, bam, bai ]
        ch_slice_inputs = WorkflowLilac.getSliceInputs(ch_lilac_bams, tumor_sequence_type)
        // Isolate meta containing expected file count to use for non-blocking groupTuple later
        ch_slice_meta_individual = ch_slice_inputs
            .map {
                def meta_lilac = it[0]
                return [key: meta_lilac.key, count: meta_lilac.count]
            }
        // Apply slicing to unique files only
        // channel: [ meta_lilac, bam, bai ]
        ch_slice_inputs_unique = WorkflowLilac.getUniqueInputFiles(ch_slice_inputs)
        SLICEBAM(
            ch_slice_inputs_unique,
            hla_slice_bed,
        )
        ch_versions = ch_versions.mix(SLICEBAM.out.versions)

        // Realign contigs if using 38 (use of ALT contigs implied)
        // channel: [ meta_lilac, bam, bai ]
        ch_slices_out = SLICEBAM.out.bam
        if (params.ref_data_genome_type == 'alt') {
            // Align reads with chr6
            // NOTE(SW): the aim of this process is to take reads mapping to ALT contigs and align them
            // to the three relevant HLA genes on chr6. All reads including those previously mapped to chr6
            // are realigned for consistency.
            EXTRACTCONTIG(
                'chr6',
                genome_fasta,
                genome_fai,
            )
            ch_versions = ch_versions.mix(EXTRACTCONTIG.out.versions)

            ch_slice_bams = SLICEBAM.out.bam
                .branch { meta_lilac, bam, bai ->
                    def sequence_type = Utils.getEnumFromString(meta_lilac.sequence_type_str, Constants.SequenceType)
                    wts: sequence_type == Constants.SequenceType.WTS
                    // WGS or targetted
                    non_wts: sequence_type != Constants.SequenceType.WTS
                }

            REALIGNREADS(
                ch_slice_bams.non_wts,
                EXTRACTCONTIG.out.contig,
                EXTRACTCONTIG.out.bwa_indices,
            )
            ch_slices_out = Channel.empty().mix(REALIGNREADS.out.bam, ch_slice_bams.wts)
            ch_versions = ch_versions.mix(REALIGNREADS.out.versions)
        }

        // Re-replicate and flow expected file count into meta
        // channel: [ meta_lilac, [sequence_type_str, sample_type_str, bam, bai] ]
        ch_slices_out_individual = ch_slices_out
            .flatMap{ data ->
                def meta_lilac = data[0]
                def fps = data[1..-1]
                meta_lilac.keys.collect { key ->
                    return [[key: key], [meta_lilac.sequence_type_str, meta_lilac.sample_type_str, *fps]]
                }
            }
        // Adding expected file count
        // channel: [ meta_lilac, [sequence_type_str, sample_type_str, bam, bai] ]
        ch_slices_ready = WorkflowOncoanalyser.joinMeta(
            ch_slices_out_individual,
            ch_slice_meta_individual,
            key_a: 'key',
        )

        // Gather and order files from same grouping using non-blocked groupTuple via provided file counts
        // channel: [ meta_lilac, normal_wgs_bam, normal_wgs_bai, tumor_bam, tumor_bai, tumor_wts_bam, tumor_wts_bai ]
        ch_slices_organised = WorkflowLilac.sortSlices(ch_slices_ready, tumor_sequence_type)

        // Restore original meta so we can join with PURPLE directory
        // channel: [val(meta)]
        ch_metas = ch_lilac_bams.map { return it[0] }
        // channel: [val(meta), normal_wgs_bam, normal_wgs_bai, tumor_bam, tumor_bai, tumor_wts_bam, tumor_wts_bai]
        ch_lilac_inputs_slices = WorkflowOncoanalyser.restoreMeta(
            ch_slices_organised,
            ch_metas,
        )

        // Get inputs from PURPLE
        // channel: [ meta, gene_cn ]
        ch_lilac_inputs_gene_cn = ch_lilac_inputs_purple
            .map { meta, purple_dir ->
                if (purple_dir == []) {
                    return [meta, []]
                }

                def tumor_id = Utils.getTumorSampleName(meta, run_config.mode)
                def gene_cn = file(purple_dir).resolve("${tumor_id}.purple.cnv.gene.tsv")
                return gene_cn.exists() ? [meta, gene_cn] : [meta, []]
            }

        // channel: [ meta, smlv_vcf ]
        ch_lilac_inputs_smlv_vcf = ch_lilac_inputs_purple
            .map { meta, purple_dir ->
                if (purple_dir == []) {
                    return [meta, []]
                }

                def tumor_id = Utils.getTumorSampleName(meta, run_config.mode)
                def smlv_vcf = file(purple_dir).resolve("${tumor_id}.purple.somatic.vcf.gz")
                return smlv_vcf.exists() ? [meta, smlv_vcf] : [meta, []]
            }

        // channel: [ meta, normal_wgs_bam, normal_wgs_bai, tumor_bam, tumor_bai, tumor_wts_bam, tumor_wts_bai, purple_dir ]
        ch_lilac_inputs_full = WorkflowOncoanalyser.groupByMeta(
            ch_lilac_inputs_slices,
            ch_lilac_inputs_gene_cn,
            ch_lilac_inputs_smlv_vcf,
            flatten_mode: 'nonrecursive',
        )

        // Create final input channel for LILAC, remove samples with only WTS BAMs
        // channel: [ meta_lilac, normal_wgs_bam, normal_wgs_bai, tumor_bam, tumor_bai, tumor_wts_bam, tumor_wts_bai, gene_cn, smlv_vcf ]
        ch_lilac_inputs = ch_lilac_inputs_full
            .map {
                def meta = it[0]
                def fps = it[1..-1]

                // LILAC requires either tumor or normal WGS BAM
                if (fps[0] == [] && fps[2] == []) {
                    return Constants.PLACEHOLDER_META
                }

                def normal_id_key = ['sample_name', Constants.SampleType.NORMAL, Constants.SequenceType.WGS]
                def tumor_id_key = ['sample_name', Constants.SampleType.TUMOR, tumor_sequence_type]

                def meta_lilac = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: meta.containsKey(tumor_id_key) ? meta.getAt(tumor_id_key) : '',
                    normal_id: meta.containsKey(normal_id_key) ? meta.getAt(normal_id_key) : '',
                ]
                return [meta_lilac, *fps]
            }
            .filter { it != Constants.PLACEHOLDER_META }

        // Run LILAC
        LILAC(
            ch_lilac_inputs,
            genome_fasta,
            genome_version,
            lilac_resource_dir,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(LILAC.out.versions)
        ch_outputs = WorkflowOncoanalyser.restoreMeta(LILAC.out.lilac_dir, ch_inputs)

    emit:
        lilac_dir = ch_outputs // channel: [ meta, lilac_dir ]

        versions = ch_versions // channel: [ versions.yml ]
}
