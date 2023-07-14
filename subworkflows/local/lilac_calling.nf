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

        // Create channels for available input BAMs
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

        } else if (run_config.mode == Constants.RunMode.PANEL) {

            ch_lilac_bams = ch_inputs
                .map { meta ->
                    def tumor_bam = Utils.getTumorBam(meta, run_config.mode)
                    [meta, tumor_bam, [], "${tumor_bam}.bai", []]
                }

        } else {
            ch_lilac_bams = ch_inputs.map { meta -> [meta, [], [], [], []] }
        }

        // Set tumor sequence type for non-WTS input (i.e. WGS or panel)
        switch(run_config.mode) {
            case Constants.RunMode.WGS:
            case Constants.RunMode.WGTS:
                tumor_sequence_type = Constants.SequenceType.WGS
                break
            case Constants.RunMode.PANEL:
                tumor_sequence_type = Constants.SequenceType.TARGETTED
                break
            default:
                assert false
        }

        // Realign reads mapping to HLA regions and homologus regions if using reference genome with ALT contigs
        // NOTE(SW): the aim of this process is to take reads mapping to ALT contigs and align them to the three
        // relevant HLA genes on chr6. All reads including those previously mapped to chr6 are realigned for
        // consistency.
        // channel: [ meta, bam, bai ]
        ch_lilac_bams = ch_lilac_bams
        if (params.ref_data_genome_type == 'alt') {

            // Split non-WTS BAMs into tumor and normal, accounting for optional input
            // channel: [ meta_extra, bam, bai ]
            ch_slice_bams = ch_lilac_bams
                .flatMap { meta, tumor_bam, normal_bam, tumor_bai, normal_bai ->
                    return [
                        [[key: meta.id, *:meta, sequence_type: 'tumor'], tumor_bam, tumor_bai],
                        [[key: meta.id, *:meta, sequence_type: 'normal'], normal_bam, normal_bai],
                    ]
                }
                .branch {
                    def empty = it[1..-1] == [[], []]
                    present: !empty
                    absent: empty
                }

            // Prepare slice input channel
            // channel: [ meta_slice, bam, bai ]
            ch_slice_bams_input = ch_slice_bams.present
                .map { meta, bam, bai ->

                    def sample_name
                    if (meta.sequence_type == 'tumor') {
                        sample_name = meta[['sample_name', Constants.SampleType.TUMOR, tumor_sequence_type]]
                    } else if (meta.sequence_type == 'normal') {
                        sample_name = meta[['sample_name', Constants.SampleType.NORMAL, Constants.SequenceType.WGS]]
                    }

                    def meta_slice = [
                        key: meta.id,
                        id: sample_name,
                        sequence_type: meta.sequence_type,
                    ]

                    return [meta_slice, bam, bai]

                }

            //
            // MODULE: Custom BAM slice (LILAC)
            //
            SLICEBAM(
                ch_slice_bams_input,
                hla_slice_bed,
            )
            ch_versions = ch_versions.mix(SLICEBAM.out.versions)

            //
            // MODULE: Custom extract contig (LILAC)
            //
            // Extract chromosome 6 from the reference and create BWA indexes
            EXTRACTCONTIG(
                'chr6',
                genome_fasta,
                genome_fai,
            )
            ch_versions = ch_versions.mix(EXTRACTCONTIG.out.versions)

            //
            // MODULE: Custom realign reads (LILAC)
            //
            // Realign selected reads relevant to HLA to chromosome 6
            REALIGNREADS(
                SLICEBAM.out.bam,
                EXTRACTCONTIG.out.contig,
                EXTRACTCONTIG.out.bwa_indices,
            )
            ch_versions = ch_versions.mix(REALIGNREADS.out.versions)

            // Split realigned BAMs by sequence type for processing below
            // channel: [ meta, bam, bai ]
            ch_slice_reunited_bams = REALIGNREADS.out.bam
                .mix (ch_slice_bams.absent)
                .branch {
                    def meta = it[0]
                    tumor: meta.sequence_type == 'tumor'
                    normal: meta.sequence_type == 'normal'
                }

            // Restore original meta for realigned BAMs then combine with absent inputs
            // channel: [ meta, tumor_bam, normal_bam, tumor_bai, normal_bai ]
            ch_lilac_bams = WorkflowOncoanalyser.groupByMeta(
                WorkflowOncoanalyser.restoreMeta(ch_slice_reunited_bams.tumor, ch_inputs),
                WorkflowOncoanalyser.restoreMeta(ch_slice_reunited_bams.normal, ch_inputs),
            )
                .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
                    return [meta, tumor_bam, normal_bam, tumor_bai, normal_bai]
                }

        }

        // Combine non-WTS and WTS BAMs, and order as required
        // channel: [ meta, normal_wgs_bam, normal_wgs_bai, tumor_bam, tumor_bai, tumor_wts_bam, tumor_wts_bai ]
        ch_lilac_bams_combined = WorkflowOncoanalyser.groupByMeta(
            ch_lilac_bams_wts,
            ch_lilac_bams,
            flatten_mode: 'nonrecursive',
        )
            .map { data ->
                def meta = data[0]
                def (tbam_wts, tbai_wts, tbam, nbam, tbai, nbai) = data[1..-1]
                return [meta, nbam, nbai, tbam, tbai, tbam_wts, tbai_wts]
            }

        // Add PURPLE inputs
        // channel: [ meta, normal_wgs_bam, normal_wgs_bai, tumor_bam, tumor_bai, tumor_wts_bam, tumor_wts_bai, purple_dir ]
        ch_lilac_inputs_full = WorkflowOncoanalyser.groupByMeta(
            ch_lilac_bams_combined,
            ch_purple,
            flatten_mode: 'nonrecursive',
        )

        // Create final input channel for LILAC, remove samples with only WTS BAMs
        // channel: [ meta_lilac, normal_wgs_bam, normal_wgs_bai, tumor_bam, tumor_bai, tumor_wts_bam, tumor_wts_bai, purple_dir ]
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
