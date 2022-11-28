//
// SV Prep is a BAM filter designed to select only reads relevant to SV events run prior to GRIDSS.
// GRIDSS is a software suite containing tools useful for the detection of genomic rearrangements.
//

include { ASSEMBLE        as GRIDSS_ASSEMBLE        } from '../../modules/local/svprep/assemble/main'
include { CALL            as GRIDSS_CALL            } from '../../modules/local/svprep/call/main'
include { PREPROCESS      as GRIDSS_PREPROCESS      } from '../../modules/local/svprep/preprocess/main'
include { SVPREP          as SVPREP_NORMAL          } from '../../modules/local/svprep/svprep/main'
include { SVPREP          as SVPREP_TUMOR           } from '../../modules/local/svprep/svprep/main'
include { DEPTH_ANNOTATOR as SVPREP_DEPTH_ANNOTATOR } from '../../modules/local/svprep/depth_annotator/main'

workflow GRIDSS_SVPREP {
  take:
    ch_inputs                       // channel: [val(meta)]
    gridss_config                   //    file: /path/to/gridss_config (optional)
    ref_data_genome_fa              //    file: /path/to/genome_fa
    ref_data_genome_version         //     val: genome version
    ref_data_genome_fai             //    file: /path/to/genome_fai
    ref_data_genome_dict            //    file: /path/to/genome_dict
    ref_data_genome_bwa_index       //    file: /path/to/genome_bwa_index_dir/
    ref_data_genome_bwa_index_image //    file: /path/to/genome_bwa_index_image
    ref_data_genome_gridss_index    //    file: /path/to/genome_gridss_index
    ref_data_gridss_blacklist       //     val: /path/to/gridss_blacklist
    ref_data_sv_prep_blacklist      //    file: /path/to/sv_prep_blacklist
    ref_data_known_fusions          //    file: /path/to/known_fusions

  main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Prepare tumor sample inputs
    // channel: [val(meta_svprep), bam, bai, []]
    ch_svprep_tumor_inputs = ch_inputs
      .map { meta ->
        def meta_svprep = [id: meta.get(['sample_name', 'tumor']), meta_full: meta]
        def bam = meta.get(['bam_wgs', 'tumor'])
        def bai = "${bam}.bai"
        return [meta_svprep, bam, bai, []]
      }

    // Filter tumor BAM
    SVPREP_TUMOR(
      ch_svprep_tumor_inputs,
      ref_data_genome_fa,
      ref_data_genome_version,
      ref_data_sv_prep_blacklist,
      ref_data_known_fusions,
      'JUNCTIONS;BAM;FRAGMENT_LENGTH_DIST', // -write_types argument switch and value
      false,                                // -calc_fragment_length argument switch
    )
    ch_versions = ch_versions.mix(SVPREP_TUMOR.out.versions)

    // Prepare normal sample inputs
    // channel: [val(meta_svprep), bam_normal, bai_normal, junctions_tumor]
    ch_svprep_normal_inputs = SVPREP_TUMOR.out.junctions
      .map { meta_tumor_svprep, junctions_tumor ->
        def meta = meta_tumor_svprep['meta_full']
        def meta_svprep = [id: meta.get(['sample_name', 'normal']), meta_full: meta]
        def bam_normal = meta.get(['bam_wgs', 'normal'])
        return [meta_svprep, bam_normal, "${bam_normal}.bai", junctions_tumor]
      }

    SVPREP_NORMAL(
      ch_svprep_normal_inputs,
      ref_data_genome_fa,
      ref_data_genome_version,
      ref_data_sv_prep_blacklist,
      ref_data_known_fusions,
      false, // -write_types argument switch and value
      true,  // -calc_fragment_length argument switch
    )
    ch_versions = ch_versions.mix(SVPREP_NORMAL.out.versions)

    ch_preprocess_inputs_tumor = WorkflowHmftools.group_by_meta(
      SVPREP_TUMOR.out.bam,
      ch_svprep_tumor_inputs,
    )
      .map { meta_svprep, bam_filtered, bam, bai ->
        def meta_gridss = [
          id: meta_svprep.id,
          sample_type: 'tumor',
          subject_id: meta_svprep.meta_full.id,
        ]
        return [meta_gridss, bam, bam_filtered]
      }

    ch_preprocess_inputs_normal = WorkflowHmftools.group_by_meta(
      SVPREP_NORMAL.out.bam,
      ch_svprep_normal_inputs,
    )
      .map { meta_svprep, bam_filtered, bam, bai, junctions ->
        def meta_gridss = [
          id: meta_svprep.id,
          sample_type: 'normal',
          subject_id: meta_svprep.meta_full.id,
        ]
        return [meta_gridss, bam, bam_filtered]
      }

    ch_preprocess_inputs = Channel.empty()
      .mix(
        ch_preprocess_inputs_tumor,
        ch_preprocess_inputs_normal,
      )

    // Preprocess reads
    GRIDSS_PREPROCESS(
      ch_preprocess_inputs,
      gridss_config,
      ref_data_genome_fa,
      ref_data_genome_fai,
      ref_data_genome_dict,
      ref_data_genome_bwa_index,
      ref_data_genome_bwa_index_image,
      ref_data_genome_gridss_index,
    )
    ch_versions = ch_versions.mix(GRIDSS_PREPROCESS.out.versions)

    // Gather BAMs and outputs from preprocessing for each tumor/normal set
    // channel: [subject_id, [[val(meta_gridss), bam, bam_filtered, preprocess_dir], ...]]
    ch_bams_and_preprocess = WorkflowHmftools.group_by_meta(
      ch_preprocess_inputs,
      GRIDSS_PREPROCESS.out.preprocess_dir,
    )
    .map { [it[0].subject_id, it] }
    .groupTuple(size: 2)

    // Order and organise inputs for assembly
    // channel: [val(meta_gridss), [bams], [bams_filtered], [preprocess_dirs], [labels]]
    ch_assemble_inputs = ch_bams_and_preprocess
      .map { subject_id, entries ->
        def (tmeta, tbam, tbam_filtered, tpreprocess) = entries.find { e -> e[0].sample_type == 'tumor' }
        def (nmeta, nbam, nbam_filtered, npreprocess) = entries.find { e -> e[0].sample_type == 'normal' }
        def meta_gridss = [id: tmeta.subject_id]
        return [
          meta_gridss,
          [nbam, tbam],
          [nbam_filtered, tbam_filtered],
          [npreprocess, tpreprocess],
          [nmeta.id, tmeta.id],
        ]
      }

    // Assemble variants
    GRIDSS_ASSEMBLE(
      ch_assemble_inputs,
      gridss_config,
      ref_data_genome_fa,
      ref_data_genome_fai,
      ref_data_genome_dict,
      ref_data_genome_bwa_index,
      ref_data_genome_bwa_index_image,
      ref_data_genome_gridss_index,
      ref_data_gridss_blacklist,
    )
    ch_versions = ch_versions.mix(GRIDSS_ASSEMBLE.out.versions)

    // Prepare inputs for variant calling
    // channel: [val(meta_gridss), [bams], [bams_filtered], assemble_dir, [labels]]
    ch_call_inputs = WorkflowHmftools.group_by_meta(
      ch_assemble_inputs,
      GRIDSS_ASSEMBLE.out.assemble_dir,
      flatten: false,
    )
      .map { data ->
        def meta = data[0]
        def (bams, bams_filtered, preprocess_dirs, labels) = data[1]
        def (assemble_dir) = data[2]
        return [meta, bams, bams_filtered, assemble_dir, labels]
      }

    // Call variants
    GRIDSS_CALL(
      ch_call_inputs,
      gridss_config,
      ref_data_genome_fa,
      ref_data_genome_fai,
      ref_data_genome_dict,
      ref_data_genome_bwa_index,
      ref_data_genome_bwa_index_image,
      ref_data_genome_gridss_index,
      ref_data_gridss_blacklist,
    )
    ch_versions = ch_versions.mix(GRIDSS_CALL.out.versions)

    // Prepare inputs for depth annotation
    // channel: [val(meta), [bams], [bais], vcf, [labels]]
    ch_depth_inputs = WorkflowHmftools.group_by_meta(
      ch_inputs.map { meta -> [meta.id, meta] },
      GRIDSS_CALL.out.vcf.map { meta, vcf -> [meta.id, vcf] },
    )
      .map { id, meta, vcf ->
        def tbam = meta.get(['bam_wgs', 'tumor'])
        def nbam = meta.get(['bam_wgs', 'normal'])
        return [
          meta,
          [nbam, tbam],
          ["${nbam}.bai", "${tbam}.bai"],
          vcf,
          [meta.get(['sample_name', 'normal']), meta.get(['sample_name', 'tumor'])],
        ]
      }

    // Add depth annotations to SVs
    SVPREP_DEPTH_ANNOTATOR(
      ch_depth_inputs,
      ref_data_genome_fa,
      ref_data_genome_version,
    )

  emit:
    results  = SVPREP_DEPTH_ANNOTATOR.out.vcf // channel: [val(meta), vcf]

    versions = ch_versions                    // channel: [versions.yml]
}
