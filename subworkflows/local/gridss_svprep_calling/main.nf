//
// SV Prep selects only reads relevant to SV events run prior to execution of GRIDSS.
// GRIDSS detects structural variants, and reports breakends and breakpoints.
//

import Constants
import Utils

include { GRIDSS_ASSEMBLE as ASSEMBLE               } from '../../../modules/local/svprep/assemble/main'
include { GRIDSS_CALL as CALL                       } from '../../../modules/local/svprep/call/main'
include { SVPREP_DEPTH_ANNOTATOR as DEPTH_ANNOTATOR } from '../../../modules/local/svprep/depth_annotator/main'
include { GRIDSS_PREPROCESS as PREPROCESS           } from '../../../modules/local/svprep/preprocess/main'
include { SVPREP as SVPREP_NORMAL                   } from '../../../modules/local/svprep/svprep/main'
include { SVPREP as SVPREP_TUMOR                    } from '../../../modules/local/svprep/svprep/main'

workflow GRIDSS_SVPREP_CALLING {
    take:
        // Sample data
        ch_inputs              // channel: [mandatory] [ meta ]
        ch_tumor_bam           // channel: [mandatory] [ meta, bam, bai ]
        ch_normal_bam          // channel: [mandatory] [ meta, bam, bai ]

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

        // Select input sources and sort
        // channel: runnable_tn: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai ]
        // channel: runnable_to: [ meta, tumor_bam, tumor_bai ]
        // channel: skip: [ meta ]
        ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
            ch_tumor_bam,
            ch_normal_bam,
        )
            .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
                return [
                    meta,
                    Utils.selectCurrentOrExisting(tumor_bam, meta, Constants.INPUT.BAM_MARKDUPS_DNA_TUMOR),
                    tumor_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_TUMOR),
                    Utils.selectCurrentOrExisting(normal_bam, meta, Constants.INPUT.BAM_MARKDUPS_DNA_NORMAL),
                    normal_bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_NORMAL),
                ]
            }
            .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->

                def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.GRIDSS_VCF)

                runnable_tn: tumor_bam && normal_bam && !has_existing
                runnable_to: tumor_bam && !has_existing
                    return [meta, tumor_bam, tumor_bai]
                skip: true
                    return meta
            }

        //
        // MODULE: SV Prep (tumor)
        //
        // Create process input channel
        // channel: [ meta_svprep, bam_tumor, bai_tumor, [] ]
        ch_svprep_tumor_inputs = Channel.empty()
            .mix(
                ch_inputs_sorted.runnable_to.map { [*it, [], []] },
                ch_inputs_sorted.runnable_tn,
            )
            .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->

                def meta_svprep = [
                    key: meta.group_id,
                    id: meta.group_id,
                    sample_id: Utils.getTumorDnaSampleName(meta),
                    sample_type: 'tumor',
                    // NOTE(SW): slightly redundant since we have this information then lose it with .mix above
                    group_size: normal_bam ? 2 : 1
                ]

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
        // channel: [ meta_svprep, bam_normal, bai_normal, junctions_tumor ]
        ch_svprep_normal_inputs = WorkflowOncoanalyser.groupByMeta(
            ch_inputs_sorted.runnable_tn,
            // NOTE(SW): this implicitly selects only entries present in ch_inputs_sorted.runnable_tn
            WorkflowOncoanalyser.restoreMeta(SVPREP_TUMOR.out.junctions, ch_inputs_sorted.runnable_tn.map { it[0] })
        )
            .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, junctions_tumor ->

                def meta_svprep = [
                    key: meta.group_id,
                    id: meta.group_id,
                    sample_id: Utils.getNormalDnaSampleName(meta),
                    sample_type: 'normal',
                    group_size: 2,  // Assumption holds since germline only is not supported and we source from runnable_tn
                ]

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
            .map { meta_svprep, bam, bam_filtered ->

                def meta_gridss = [
                    key: meta_svprep.key,
                    id: "${meta_svprep.id}__${meta_svprep.sample_id}",
                    sample_id: meta_svprep.sample_id,
                    sample_type: meta_svprep.sample_type,
                    group_size: meta_svprep.group_size,
                ]

                return [meta_gridss, bam, bam_filtered]
            }

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
                        tmeta.sample_id,
                    ]

                } else {

                    def (nmeta, nbam, nbam_filtered, npreprocess) = normal_entry
                    data = [
                        meta_gridss,
                        [nbam, tbam],
                        [nbam_filtered, tbam_filtered],
                        [npreprocess, tpreprocess],
                        [nmeta.sample_id, tmeta.sample_id],
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
        // channel: [ meta, [bams], [bais], vcf, [labels] ]
        ch_depth_inputs_tn = WorkflowOncoanalyser.groupByMeta(
            ch_inputs_sorted.runnable_tn,
            // NOTE(SW): this implicitly selects only entries present in ch_inputs_sorted.runnable_tn
            WorkflowOncoanalyser.restoreMeta(CALL.out.vcf, ch_inputs_sorted.runnable_tn.map { it[0] })
        )
            .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, vcf  ->
                return [
                    meta,
                    [normal_bam, tumor_bam],
                    [normal_bai, tumor_bai],
                    vcf,
                    [Utils.getNormalDnaSampleName(meta), Utils.getTumorDnaSampleName(meta)],
                ]
            }

        // channel: [ meta, bam, bai, vcf, label ]
        ch_depth_inputs_to = WorkflowOncoanalyser.groupByMeta(
            ch_inputs_sorted.runnable_to,
            // NOTE(SW): this implicitly selects only entries present in ch_inputs_sorted.runnable_to
            WorkflowOncoanalyser.restoreMeta(CALL.out.vcf, ch_inputs_sorted.runnable_to.map { it[0] })
        )
            .map { meta, tumor_bam, tumor_bai, vcf ->
                return [
                    meta,
                    tumor_bam,
                    tumor_bai,
                    vcf,
                    Utils.getTumorDnaSampleName(meta),
                ]
            }

        // channel: runnable_tn: [ meta_svprep, [bams], [bais], vcf, [labels] ]
        // channel: runnable_to: [ meta_svprep, bam, bai, vcf, label ]
        ch_depth_inputs = Channel.empty()
            .mix(
                ch_depth_inputs_tn,
                ch_depth_inputs_to,
            )
            .map { d ->

                def meta = d[0]
                def fps = d[1..-1]

                def meta_svprep = [
                    key: meta.group_id,
                    id: meta.group_id,
                    tumor_id: Utils.getTumorDnaSampleName(meta)
                ]

                return [meta_svprep, *fps]
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
