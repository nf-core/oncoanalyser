//
// PURPLE is a CNV caller that infers purity/ploidy and recovers low-confidence SVs
//

include { PURPLE } from '../../modules/local/purple/main'

workflow PURPLE_CALLING {
    take:
        // Sample data
        ch_inputs
        ch_amber
        ch_cobalt
        ch_smlv_somatic
        ch_smlv_germline
        ch_sv_somatic
        ch_sv_somatic_recovery

        // Reference data
        ref_data_genome_fasta
        ref_data_genome_fai
        ref_data_genome_dict
        ref_data_genome_version
        gc_profile
        sage_known_hotspots_somatic
        sage_known_hotspots_germline
        driver_gene_panel
        ensembl_data_resources
        purple_germline_del

        // Params
        run

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [val(meta), sv_vcf, sv_tbi, sv_recovery_vcf, sv_recovery_tbi]
        ch_purple_inputs_sv = Channel.empty()
        if (run.gripss) {
            // NOTE(SW): GRIPSS will be run for all WGS entries, so no optionals here
            ch_purple_inputs_sv = WorkflowOncoanalyser.groupByMeta(
                ch_sv_somatic,
                ch_sv_somatic_recovery,
            )
        } else {
            ch_purple_inputs_sv = WorkflowOncoanalyser.groupByMeta(
                WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.GRIPSS_VCF_TUMOR, type: 'optional'),
                WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.GRIPSS_UNFILTERED_VCF_TUMOR, type: 'optional'),
                flatten_mode: 'nonrecursive',
            )
              .map { meta, vcf, vcf_unfiltered->
                  def tbi            = vcf == []            ? [] : "${vcf}.tbi"
                  def tbi_unfiltered = vcf_unfiltered == [] ? [] : "${vcf_unfiltered}.tbi"
                  return [meta, vcf, tbi, vcf_unfiltered, tbi_unfiltered]
              }
        }

        // channel: [val(meta), amber_dir, cobalt_dir, sv_vcf, sv_tbi, sv_recovery_vcf, sv_recovery_tbi, smlv_tumor_vcf, smlv_normal_vcf]
        ch_purple_inputs_source = WorkflowOncoanalyser.groupByMeta(
            // Required inputs
            run.amber ? ch_amber: WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.AMBER_DIR),
            run.cobalt ? ch_cobalt: WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.COBALT_DIR),
            // Optional inputs
            ch_purple_inputs_sv,
            run.pave ? ch_smlv_somatic : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PAVE_VCF_TUMOR, type: 'optional'),
            run.pave ? ch_smlv_germline : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PAVE_VCF_NORMAL, type: 'optional'),
            flatten_mode: 'nonrecursive',
        )

        // Create process-specific meta
        // channel: [val(meta_purple), amber_dir, cobalt_dir, sv_vcf, sv_tbi, sv_recovery_vcf, sv_recovery_tbi, smlv_tumor_vcf, smlv_normal_vcf]
        ch_purple_inputs = ch_purple_inputs_source
            .map {
                def meta = it[0]
                def other = it[1..-1]
                def meta_purple = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: Utils.getTumorWgsSampleName(meta),
                    normal_id: Utils.getNormalWgsSampleName(meta),
              ]
              return [meta_purple, *other]
            }

        PURPLE(
            ch_purple_inputs,
            ref_data_genome_fasta,
            ref_data_genome_fai,
            ref_data_genome_dict,
            ref_data_genome_version,
            gc_profile,
            sage_known_hotspots_somatic,
            sage_known_hotspots_germline,
            driver_gene_panel,
            ensembl_data_resources,
            purple_germline_del,
        )

        ch_outputs = WorkflowOncoanalyser.restoreMeta(PURPLE.out.purple_dir, ch_inputs)
        ch_versions = ch_versions.mix(PURPLE.out.versions)

    emit:
        purple_dir = ch_outputs  // channel: [val(meta), purple_dir]

        versions   = ch_versions // channel: [versions.yml]
}
