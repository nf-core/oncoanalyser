//
// SV Prep selects only reads relevant to SV events run prior to execution of GRIDSS.
// GRIDSS detects structural variants, and reports breakends and breakpoints.
//

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
        run_config             // channel: [mandatory] run configuration

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        //
        // MODULE: SV Prep (tumor)
        //
        // Prepare tumor sample inputs
        // channel: [ meta_svprep, bam_tumor, bai_tumor, [] ]
        ch_svprep_tumor_inputs = ch_inputs
            .map { meta ->
                def meta_svprep = [
                    key: meta.id,
                    id: Utils.getTumorSampleName(meta, run_config.mode),
                    sample_type: 'tumor',
                ]
                def tumor_bam = Utils.getTumorBam(meta, run_config.mode)
                return [meta_svprep, tumor_bam, "${tumor_bam}.bai", []]
            }

        // Filter tumor BAM
        SVPREP_TUMOR(
            ch_svprep_tumor_inputs,
            genome_fasta,
            genome_version,
            sv_prep_blocklist,
            known_fusions,
            'JUNCTIONS;BAM;FRAGMENT_LENGTH_DIST', // -write_types argument switch and value
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
        // channel: [ meta_gridss, bam_normal, bam_normal_filtered ]
        ch_preprocess_inputs_normal = Channel.empty()
        if (run_config.type == Constants.RunType.TUMOR_NORMAL) {

            assert [Constants.RunMode.WGS, Constants.RunMode.WGTS].contains(run_config.mode)

            // Prepare normal sample inputs
            // channel: [val(meta_svprep), bam_normal, bai_normal, junctions_tumor]
            ch_svprep_normal_inputs = WorkflowOncoanalyser.restoreMeta(SVPREP_TUMOR.out.junctions, ch_inputs)
                .map { meta, junctions_tumor ->
                    def meta_svprep = [
                        key: meta.id,
                        id: Utils.getNormalWgsSampleName(meta),
                        sample_type: 'normal',
                    ]
                    def normal_bam = Utils.getNormalWgsBam(meta)
                    return [meta_svprep, normal_bam, "${normal_bam}.bai", junctions_tumor]
                }

            SVPREP_NORMAL(
                ch_svprep_normal_inputs,
                genome_fasta,
                genome_version,
                sv_prep_blocklist,
                known_fusions,
                false, // -write_types argument switch and value
            )
            ch_versions = ch_versions.mix(SVPREP_NORMAL.out.versions)

            // channel: [ meta_gridss, bam_normal, bam_normal_filtered ]
            ch_preprocess_inputs_normal = WorkflowOncoanalyser.groupByMeta(
                SVPREP_NORMAL.out.bam,
                ch_svprep_normal_inputs,
            )
                .map { meta_svprep, bam_filtered, bam, bai, junctions ->
                    return [meta_svprep, bam, bam_filtered]
                }
        }

        //
        // MODULE: GRIDSS preprocess
        //
        // channel: [ meta_gridss, bam, bam_filtered ]
        ch_preprocess_inputs = Channel.empty()
            .mix(
                ch_preprocess_inputs_tumor,
                ch_preprocess_inputs_normal,
            )

        // Preprocess reads
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

        // Gather BAMs and outputs from preprocessing for each tumor/normal set
        // channel: [key, [[ meta_gridss, bam, bam_filtered, preprocess_dir], ...] ]
        preprocess_group_tuple_size = run_config.type == Constants.RunType.TUMOR_NORMAL ? 2 : 1
        ch_bams_and_preprocess = WorkflowOncoanalyser.groupByMeta(
            ch_preprocess_inputs,
            PREPROCESS.out.preprocess_dir,
        )
        .map { [it[0].key, it] }
        .groupTuple(size: preprocess_group_tuple_size)

        //
        // MODULE: GRIDSS assemble
        //
        // Order and organise inputs for assembly
        // channel: [ meta_gridss, [bams], [bams_filtered], [preprocess_dirs], [labels] ]
        ch_assemble_inputs = ch_bams_and_preprocess
            .map { key, entries ->
                def (tmeta, tbam, tbam_filtered, tpreprocess) = entries.find { e -> e[0].sample_type == 'tumor' }
                def meta_gridss = [id: tmeta.key]

                def data = []
                if (run_config.type == Constants.RunType.TUMOR_ONLY) {

                    data = [
                        meta_gridss,
                        tbam,
                        tbam_filtered,
                        tpreprocess,
                        tmeta.id,
                    ]

                } else if (run_config.type == Constants.RunType.TUMOR_NORMAL) {

                    def (nmeta, nbam, nbam_filtered, npreprocess) = entries.find { e -> e[0].sample_type == 'normal' }
                    data = [
                        meta_gridss,
                        [nbam, tbam],
                        [nbam_filtered, tbam_filtered],
                        [npreprocess, tpreprocess],
                        [nmeta.id, tmeta.id],
                    ]

                } else {
                    assert false
                }

                return data
            }

        // Assemble variants
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
        // Prepare inputs for variant calling
        // channel: [ meta_gridss, [bams], [bams_filtered], assemble_dir, [labels] ]
        ch_call_inputs = WorkflowOncoanalyser.groupByMeta(
            ch_assemble_inputs,
            ASSEMBLE.out.assemble_dir,
            flatten: false,
        )
            .map { data ->
                def meta = data[0]
                def (bams, bams_filtered, preprocess_dirs, labels) = data[1]
                def (assemble_dir) = data[2]
                return [meta, bams, bams_filtered, assemble_dir, labels]
            }

        // Call variants
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
        // Prepare inputs for depth annotation, restore original meta
        // channel: [ meta_svprep, [bams], [bais], vcf, [labels] ]
        ch_depth_inputs = WorkflowOncoanalyser.groupByMeta(
            ch_inputs.map { meta -> [meta.id, meta] },
            CALL.out.vcf.map { meta, vcf -> [meta.id, vcf] },
        )
            .map { id, meta, vcf ->
                def tbam = Utils.getTumorBam(meta, run_config.mode)
                def meta_svprep = [
                    id: meta.id,
                    tumor_id: Utils.getTumorSampleName(meta, run_config.mode),
                ]

                def data = []
                if (run_config.type == Constants.RunType.TUMOR_ONLY) {

                    data = [
                        meta_svprep,
                        tbam,
                        "${tbam}.bai",
                        vcf,
                        Utils.getTumorSampleName(meta, run_config.mode),
                    ]

                } else if (run_config.type == Constants.RunType.TUMOR_NORMAL) {

                    def nbam = Utils.getNormalWgsBam(meta)

                    data = [
                        meta_svprep,
                        [nbam, tbam],
                        ["${nbam}.bai", "${tbam}.bai"],
                        vcf,
                        [Utils.getNormalWgsSampleName(meta), Utils.getTumorWgsSampleName(meta)],
                    ]

                } else {
                    assert false
                }

                return data

            }

        // Add depth annotations to SVs
        DEPTH_ANNOTATOR(
            ch_depth_inputs,
            genome_fasta,
            genome_version,
        )

        // Reunite final VCF with the corresponding input meta object
        ch_out = Channel.empty()
            .concat(
                ch_inputs.map { meta -> [meta.id, meta] },
                DEPTH_ANNOTATOR.out.vcf.map { meta, vcf -> [meta.id, vcf] },
            )
            .groupTuple(size: 2)
            .map { id, other -> other.flatten() }

    emit:
        results  = ch_out      // channel: [ meta, vcf ]

        versions = ch_versions // channel: [ versions.yml ]
}
