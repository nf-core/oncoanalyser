//
// SV Prep is a BAM filter designed to select only reads relevant to SV events run prior to GRIDSS.
// GRIDSS is a software suite containing tools useful for the detection of genomic rearrangements.
//
import Utils

include { GRIDSS_ASSEMBLE as ASSEMBLE               } from '../../modules/local/svprep/assemble/main'
include { GRIDSS_CALL as CALL                       } from '../../modules/local/svprep/call/main'
include { SVPREP_DEPTH_ANNOTATOR as DEPTH_ANNOTATOR } from '../../modules/local/svprep/depth_annotator/main'
include { GRIDSS_PREPROCESS as PREPROCESS           } from '../../modules/local/svprep/preprocess/main'
include { SVPREP as SVPREP_NORMAL                   } from '../../modules/local/svprep/svprep/main'
include { SVPREP as SVPREP_TUMOR                    } from '../../modules/local/svprep/svprep/main'

include { CHANNEL_GROUP_INPUTS } from './channel_group_inputs'

workflow GRIDSS_SVPREP_CALLING {
    take:
        ch_inputs                       // channel: [val(meta)]
        gridss_config                   //    file: /path/to/gridss_config (optional)
        ref_data_genome_fasta           //    file: /path/to/genome_fasta
        ref_data_genome_version         //     val: genome version
        ref_data_genome_fai             //    file: /path/to/genome_fai
        ref_data_genome_dict            //    file: /path/to/genome_dict
        ref_data_genome_bwa_index       //    file: /path/to/genome_bwa_index_dir/
        ref_data_genome_bwa_index_image //    file: /path/to/genome_bwa_index_image
        ref_data_genome_gridss_index    //    file: /path/to/genome_gridss_index
        ref_data_gridss_blocklist       //     val: /path/to/gridss_blocklist
        ref_data_sv_prep_blocklist      //    file: /path/to/sv_prep_blocklist
        ref_data_known_fusions          //    file: /path/to/known_fusions

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Get input meta groups
        CHANNEL_GROUP_INPUTS(
            ch_inputs,
        )

        // Prepare tumor sample inputs
        // channel: [val(meta_svprep), bam_tumor, bai_tumor, []]
        ch_svprep_tumor_inputs = CHANNEL_GROUP_INPUTS.out.wgs_present
            .map { meta ->
                def meta_svprep = [
                    key: meta.id,
                    id: Utils.getTumorWgsSampleName(meta),
                    sample_type: 'tumor',
                ]
                def tumor_bam = Utils.getTumorWgsBam(meta)
                return [meta_svprep, tumor_bam, "${tumor_bam}.bai", []]
            }

        // Filter tumor BAM
        SVPREP_TUMOR(
            ch_svprep_tumor_inputs,
            ref_data_genome_fasta,
            ref_data_genome_version,
            ref_data_sv_prep_blocklist,
            ref_data_known_fusions,
            'JUNCTIONS;BAM;FRAGMENT_LENGTH_DIST', // -write_types argument switch and value
        )
        ch_versions = ch_versions.mix(SVPREP_TUMOR.out.versions)

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
            ref_data_genome_fasta,
            ref_data_genome_version,
            ref_data_sv_prep_blocklist,
            ref_data_known_fusions,
            false, // -write_types argument switch and value
        )
        ch_versions = ch_versions.mix(SVPREP_NORMAL.out.versions)

        // channel: [val(meta_gridss), bam_tumor, bam_tumor_filtered]
        ch_preprocess_inputs_tumor = WorkflowOncoanalyser.groupByMeta(
            SVPREP_TUMOR.out.bam,
            ch_svprep_tumor_inputs,
        )
            .map { meta_svprep, bam_filtered, bam, bai ->
                return [meta_svprep, bam, bam_filtered]
            }

        // channel: [val(meta_gridss), bam_normal, bam_normal_filtered]
        ch_preprocess_inputs_normal = WorkflowOncoanalyser.groupByMeta(
            SVPREP_NORMAL.out.bam,
            ch_svprep_normal_inputs,
        )
            .map { meta_svprep, bam_filtered, bam, bai, junctions ->
                return [meta_svprep, bam, bam_filtered]
            }

        // channel: [val(meta_gridss), bam, bam_filtered]
        ch_preprocess_inputs = Channel.empty()
            .mix(
                ch_preprocess_inputs_tumor,
                ch_preprocess_inputs_normal,
            )

        // Preprocess reads
        PREPROCESS(
            ch_preprocess_inputs,
            gridss_config,
            ref_data_genome_fasta,
            ref_data_genome_fai,
            ref_data_genome_dict,
            ref_data_genome_bwa_index,
            ref_data_genome_bwa_index_image,
            ref_data_genome_gridss_index,
        )
        ch_versions = ch_versions.mix(PREPROCESS.out.versions)

        // Gather BAMs and outputs from preprocessing for each tumor/normal set
        // channel: [key, [[val(meta_gridss), bam, bam_filtered, preprocess_dir], ...]]
        ch_bams_and_preprocess = WorkflowOncoanalyser.groupByMeta(
            ch_preprocess_inputs,
            PREPROCESS.out.preprocess_dir,
        )
        .map { [it[0].key, it] }
        .groupTuple(size: 2)

        // Order and organise inputs for assembly
        // channel: [val(meta_gridss), [bams], [bams_filtered], [preprocess_dirs], [labels]]
        ch_assemble_inputs = ch_bams_and_preprocess
            .map { key, entries ->
                def (tmeta, tbam, tbam_filtered, tpreprocess) = entries.find { e -> e[0].sample_type == 'tumor' }
                def (nmeta, nbam, nbam_filtered, npreprocess) = entries.find { e -> e[0].sample_type == 'normal' }
                def meta_gridss = [id: tmeta.key]
                return [
                    meta_gridss,
                    [nbam, tbam],
                    [nbam_filtered, tbam_filtered],
                    [npreprocess, tpreprocess],
                    [nmeta.id, tmeta.id],
                ]
            }

        // Assemble variants
        ASSEMBLE(
            ch_assemble_inputs,
            gridss_config,
            ref_data_genome_fasta,
            ref_data_genome_fai,
            ref_data_genome_dict,
            ref_data_genome_bwa_index,
            ref_data_genome_bwa_index_image,
            ref_data_genome_gridss_index,
            ref_data_gridss_blocklist,
        )
        ch_versions = ch_versions.mix(ASSEMBLE.out.versions)

        // Prepare inputs for variant calling
        // channel: [val(meta_gridss), [bams], [bams_filtered], assemble_dir, [labels]]
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
            ref_data_genome_fasta,
            ref_data_genome_fai,
            ref_data_genome_dict,
            ref_data_genome_bwa_index,
            ref_data_genome_bwa_index_image,
            ref_data_genome_gridss_index,
            ref_data_gridss_blocklist,
        )
        ch_versions = ch_versions.mix(CALL.out.versions)

        // Prepare inputs for depth annotation, restore original meta
        // channel: [val(meta_svprep), [bams], [bais], vcf, [labels]]
        ch_depth_inputs = WorkflowOncoanalyser.groupByMeta(
            CHANNEL_GROUP_INPUTS.out.wgs_present.map { meta -> [meta.id, meta] },
            CALL.out.vcf.map { meta, vcf -> [meta.id, vcf] },
        )
            .map { id, meta, vcf ->
                def tbam = Utils.getTumorWgsBam(meta)
                def nbam = Utils.getNormalWgsBam(meta)
                def meta_svprep = [id: meta.id]
                return [
                    meta_svprep,
                    [nbam, tbam],
                    ["${nbam}.bai", "${tbam}.bai"],
                    vcf,
                    [Utils.getNormalWgsSampleName(meta), Utils.getTumorWgsSampleName(meta)],
                ]
            }

        // Add depth annotations to SVs
        DEPTH_ANNOTATOR(
            ch_depth_inputs,
            ref_data_genome_fasta,
            ref_data_genome_version,
        )

        // Reunite final VCF with the corresponding input meta object
        ch_out = Channel.empty()
            .concat(
                CHANNEL_GROUP_INPUTS.out.wgs_present.map { meta -> [meta.id, meta] },
                DEPTH_ANNOTATOR.out.vcf.map { meta, vcf -> [meta.id, vcf] },
            )
            .groupTuple(size: 2)
            .map { id, other -> other.flatten() }

    emit:
        results  = ch_out      // channel: [val(meta), vcf]

        versions = ch_versions // channel: [versions.yml]
}
