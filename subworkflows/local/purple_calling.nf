//
// PURPLE is a CNV caller that infers purity/ploidy and recovers low-confidence SVs
//
import Constants

include { PURPLE } from '../../modules/local/purple/main'

include { CHANNEL_INPUTS_PURPLE } from './channel_inputs_purple'

workflow PURPLE_CALLING {
    take:
        // Sample data
        ch_inputs
        ch_amber
        ch_cobalt
        ch_smlv_somatic
        ch_smlv_germline
        ch_sv_somatic
        ch_sv_germline
        ch_sv_somatic_unfiltered

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
        target_region_bed
        target_region_ratios
        target_region_msi_indels

        // Params
        run_config

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [val(meta), amber_dir, cobalt_dir, sv_tumor_vcf, sv_tumor_tbi, sv_tumor_unfiltered_vcf, sv_tumor_unfiltered_tbi, smlv_tumor_vcf, smlv_normal_vcf]
        CHANNEL_INPUTS_PURPLE(
            ch_inputs,
            ch_amber,
            ch_cobalt,
            ch_smlv_somatic,
            ch_smlv_germline,
            ch_sv_somatic,
            ch_sv_germline,
            ch_sv_somatic_unfiltered,
            run_config,
        )

        // Create process-specific meta
        // channel: [val(meta_purple), amber_dir, cobalt_dir, sv_tumor_vcf, sv_tumor_tbi, sv_tumor_unfiltered_vcf, sv_tumor_unfiltered_tbi, smlv_tumor_vcf, smlv_normal_vcf]
        ch_purple_inputs = CHANNEL_INPUTS_PURPLE.out
            .map {
                def meta = it[0]
                def meta_purple = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: Utils.getTumorSampleName(meta, run_config.mode),
                ]

                if (run_config.type == Constants.RunType.TUMOR_NORMAL) {
                    meta_purple.normal_id = Utils.getNormalWgsSampleName(meta)
                }

                return [meta_purple, *it[1..-1]]
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
            target_region_bed,
            target_region_ratios,
            target_region_msi_indels,
        )

        ch_outputs = WorkflowOncoanalyser.restoreMeta(PURPLE.out.purple_dir, ch_inputs)
        ch_versions = ch_versions.mix(PURPLE.out.versions)

    emit:
        purple_dir = ch_outputs  // channel: [val(meta), purple_dir]

        versions   = ch_versions // channel: [versions.yml]
}
