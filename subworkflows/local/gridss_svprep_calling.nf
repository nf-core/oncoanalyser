//
// SV Prep selects only reads relevant to SV events run prior to execution of GRIDSS.
// GRIDSS detects structural variants, and reports breakends and breakpoints.
//

import Constants
import Utils

include { GRIDSS_ASSEMBLE as ASSEMBLE               } from '../../modules/local/svprep/assemble/main'
include { GRIDSS_CALL as CALL                       } from '../../modules/local/svprep/call/main'
include { SVPREP_DEPTH_ANNOTATOR as DEPTH_ANNOTATOR } from '../../modules/local/svprep/depth_annotator/main'
include { GRIDSS_PREPROCESS as PREPROCESS           } from '../../modules/local/svprep/preprocess/main'
include { SVPREP as SVPREP_NORMAL                   } from '../../modules/local/svprep/svprep/main'
include { SVPREP as SVPREP_TUMOR                    } from '../../modules/local/svprep/svprep/main'

workflow GRIDSS_SVPREP_CALLING {
    take:
        // Sample data
        ch_inputs              // channel: [mandatory] [ meta ]

        // Reference data
        genome_fasta           // channel: [mandatory] /path/to/genome_fasta
        genome_version         // channel: [mandatory] genome version
        genome_fai             // channel: [mandatory] /path/to/genome_fai
        genome_dict            // channel: [mandatory] /path/to/genome_dict
        genome_bwa_index       // channel: [mandatory] /path/to/genome_bwa_index_dir/
        genome_bwa_index_image // channel: [mandatory] /path/to/genome_bwa_index_image
        genome_gridss_index    // channel: [mandatory] /path/to/genome_gridss_index
        gridss_blocklist       // channel: [mandatory] /path/to/gridss_blocklist
        sv_prep_blocklist      // channel: [mandatory] /path/to/sv_prep_blocklist
        known_fusions          // channel: [mandatory] /path/to/known_fusions

        // Params
        gridss_config          // channel: [optional] /path/to/gridss_config

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Sort inputs
        ch_inputs_sorted = ch_inputs
            .branch { meta ->
                def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.GRIDSS_VCF)
                runnable_tn: Utils.hasTumorDnaBam(meta) && Utils.hasNormalDnaBam(meta) && !has_existing
                runnable_to: Utils.hasTumorDnaBam(meta) && !has_existing
                skip: true
            }

        //
        // MODULE: SV Prep (tumor)
        //
        // Create process input channel
        // channel: [ meta_svprep, bam_tumor, bai_tumor, [] ]
        ch_svprep_tumor_inputs = Channel.empty()
            .mix(
                ch_inputs_sorted.runnable_to,
                ch_inputs_sorted.runnable_tn,
            )
            .map { meta ->

                def tumor_id = Utils.getTumorDnaSampleName(meta)
                def meta_svprep = [
                    key: meta.group_id,
                    id: "${meta.group_id}__${tumor_id}",
                    sample_type: 'tumor',
                    // NOTE(SW): slightly redundant since we have this information then lose it with .mix above
                    group_size: Utils.hasNormalDnaBam(meta) ? 2 : 1
                ]

                def tumor_bam = Utils.getTumorDnaBam(meta)
                def tumor_bai = Utils.getTumorDnaBai(meta)

                return [meta_svprep, tumor_bam, tumor_bai, []]

            }

        // Run process
        SVPREP_TUMOR(
            ch_svprep_tumor_inputs,
            genome_fasta,
            genome_version,
            sv_prep_blocklist,
            known_fusions,
            'JUNCTIONS;BAM;FRAGMENT_LENGTH_DIST',  // -write_types argument
        )

        ch_versions = ch_versions.mix(SVPREP_TUMOR.out.versions)

        // channel: [ meta_gridss, bam_tumor, bam_tumor_filtered ]
        ch_preprocess_inputs_tumor = WorkflowOncoanalyser.groupByMeta(
            SVPREP_TUMOR.out.bam,
            ch_svprep_tumor_inputs,
        )
            .map { meta_svprep, bam_filtered, bam, bai, jnc_optional ->
                return [meta_svprep, bam, bam_filtered]
            }

        //
        // MODULE: SV Prep (normal)
        //
        // Create process input channel
        // NOTE(SW): this implicitly selects only entries present in ch_inputs_sorted.runnable_tn
        // channel: [ meta_svprep, bam_normal, bai_normal, junctions_tumor ]
        ch_svprep_normal_inputs = WorkflowOncoanalyser.restoreMeta(SVPREP_TUMOR.out.junctions, ch_inputs_sorted.runnable_tn)
            .map { meta, junctions_tumor ->

                def normal_id = Utils.getNormalDnaSampleName(meta)
                def meta_svprep = [
                    key: meta.group_id,
                    id: "${meta.group_id}__${normal_id}",
                    sample_type: 'normal',
                    group_size: 2,  // Assumption holds since germline only is not supported and we source from runnable_tn
                ]

                def normal_bam = Utils.getNormalDnaBam(meta)
                def normal_bai = Utils.getNormalDnaBai(meta)

                return [meta_svprep, normal_bam, normal_bai, junctions_tumor]

            }

        // Run process
        SVPREP_NORMAL(
            ch_svprep_normal_inputs,
            genome_fasta,
            genome_version,
            sv_prep_blocklist,
            known_fusions,
            'JUNCTIONS;BAM;FRAGMENT_LENGTH_DIST', // -write_types argument
        )

        ch_versions = ch_versions.mix(SVPREP_NORMAL.out.versions)

        // channel: [ meta_gridss, bam_normal, bam_normal_filtered ]
        ch_preprocess_inputs_normal = WorkflowOncoanalyser.groupByMeta(
            SVPREP_NORMAL.out.bam,
            ch_svprep_normal_inputs,
        )
            // Switching meta name here from meta_svprep
            .map { meta_gridss, bam_filtered, bam, bai, junctions ->
                return [meta_gridss, bam, bam_filtered]
            }

        //
        // MODULE: GRIDSS preprocess
        //
        // Create process input channel
        // channel: [ meta_gridss, bam, bam_filtered ]
        ch_preprocess_inputs = Channel.empty()
            .mix(
                ch_preprocess_inputs_tumor,
                ch_preprocess_inputs_normal,
            )

        // Run process
        PREPROCESS(
            ch_preprocess_inputs,
            gridss_config,
            genome_fasta,
            genome_fai,
            genome_dict,
            genome_bwa_index,
            genome_bwa_index_image,
            genome_gridss_index,
        )

        ch_versions = ch_versions.mix(PREPROCESS.out.versions)

        // Gather BAMs and outputs from preprocessing for each tumor/normal and tumor only set
        // channel: [key, [[meta_gridss, bam, bam_filtered, preprocess_dir], ...] ]
        ch_bams_and_preprocess = WorkflowOncoanalyser.groupByMeta(
            ch_preprocess_inputs,
            PREPROCESS.out.preprocess_dir,
        )
            .map {
                def meta_gridss = it[0]
                def other = it[1..-1]
                [groupKey(meta_gridss.key, meta_gridss.group_size), [meta_gridss, *other]]
            }
            .groupTuple()

        //
        // MODULE: GRIDSS assemble
        //
        // Create process input channel
        // channel: tumor/normal: [ meta_gridss, [bams], [bams_filtered], [preprocess_dirs], [labels] ]
        // channel: tumor only:   [ meta_gridss, bam, bam_filtered, preprocess_dir, label ]
        ch_assemble_inputs = ch_bams_and_preprocess
            .map { key, entries ->

                assert entries.size() == 1 || entries.size() == 2

                def tumor_entry = entries.find { e -> e[0].sample_type == 'tumor' }
                def normal_entry = entries.find { e -> e[0].sample_type == 'normal' }

                assert tumor_entry !== null

                def (tmeta, tbam, tbam_filtered, tpreprocess) = tumor_entry
                def meta_gridss = [
                    // Effectively meta.group_id, and both are required. Reminder:
                    //   * key: channel element grouping
                    //   * id: task tag
                    key: tmeta.key,
                    id: tmeta.key,
                ]

                def data = []

                if (normal_entry === null) {

                    data = [
                        meta_gridss,
                        tbam,
                        tbam_filtered,
                        tpreprocess,
                        tmeta.id,
                    ]

                } else {

                    def (nmeta, nbam, nbam_filtered, npreprocess) = normal_entry
                    data = [
                        meta_gridss,
                        [nbam, tbam],
                        [nbam_filtered, tbam_filtered],
                        [npreprocess, tpreprocess],
                        [nmeta.id, tmeta.id],
                    ]

                }

                return data
            }

        // Run process
        ASSEMBLE(
            ch_assemble_inputs,
            gridss_config,
            genome_fasta,
            genome_fai,
            genome_dict,
            genome_bwa_index,
            genome_bwa_index_image,
            genome_gridss_index,
            gridss_blocklist,
        )

        ch_versions = ch_versions.mix(ASSEMBLE.out.versions)

        //
        // MODULE: GRIDSS call
        //
        // Create process input channel
        // channel: [ meta_gridss, [bams], [bams_filtered], assemble_dir, [labels] ]
        ch_call_inputs = WorkflowOncoanalyser.groupByMeta(
            ch_assemble_inputs,
            ASSEMBLE.out.assemble_dir,
            flatten: false,
        )
            .map { data ->
                def meta_gridss = data[0]
                def (bams, bams_filtered, preprocess_dirs, labels) = data[1]
                def (assemble_dir) = data[2]
                return [meta_gridss, bams, bams_filtered, assemble_dir, labels]
            }

        // Run process
        CALL(
            ch_call_inputs,
            gridss_config,
            genome_fasta,
            genome_fai,
            genome_dict,
            genome_bwa_index,
            genome_bwa_index_image,
            genome_gridss_index,
            gridss_blocklist,
        )

        ch_versions = ch_versions.mix(CALL.out.versions)

        //
        // MODULE: SV Prep depth annotation
        //
        // Restore original meta, create process input channel
        // channel: tumor/normal: [ meta_svprep, [bams], [bais], vcf, [labels] ]
        // channel: tumor only:   [ meta_svprep, bam, bai, vcf, label ]
        ch_depth_inputs = WorkflowOncoanalyser.restoreMeta(CALL.out.vcf, ch_inputs)
            .map { meta, vcf ->

                // NOTE(SW): germline only is not currently supported
                assert Utils.hasTumorDnaBam(meta)

                def meta_svprep = [
                    // Both are required. Reminder:
                    //   * key: channel element grouping
                    //   * id: task tag
                    key: meta.group_id,
                    id: meta.group_id,
                    tumor_id: Utils.getTumorDnaSampleName(meta)
                ]

                def data = []

                if (Utils.hasNormalDnaBam(meta)) {

                    data = [
                        meta_svprep,
                        [Utils.getNormalDnaBam(meta), Utils.getTumorDnaBam(meta)],
                        [Utils.getNormalDnaBai(meta), Utils.getTumorDnaBai(meta)],
                        vcf,
                        [Utils.getNormalDnaSampleName(meta), Utils.getTumorDnaSampleName(meta)],
                    ]

                } else if (Utils.hasTumorDnaBam(meta)) {

                    data = [
                        meta_svprep,
                        Utils.getTumorDnaBam(meta),
                        Utils.getTumorDnaBai(meta),
                        vcf,
                        Utils.getTumorDnaSampleName(meta),
                    ]

                } else {
                    assert false
                }

                return data
            }

        // Add depth annotations to calls
        DEPTH_ANNOTATOR(
            ch_depth_inputs,
            genome_fasta,
            genome_version,
        )

        ch_versions = ch_versions.mix(DEPTH_ANNOTATOR.out.versions)

        // Set outputs, restoring original meta
        // channel: [ meta, gridss_vcf ]
        ch_outputs = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(DEPTH_ANNOTATOR.out.vcf, ch_inputs),
                ch_inputs_sorted.skip.map { meta -> [meta, []] },
            )

    emit:
        vcf      = ch_outputs  // channel: [ meta, vcf ]

        versions = ch_versions // channel: [ versions.yml ]
}
