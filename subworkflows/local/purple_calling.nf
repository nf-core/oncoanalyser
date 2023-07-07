//
// PURPLE is a CNV caller that infers purity/ploidy and recovers low-confidence SVs
//

import Constants

include { PURPLE } from '../../modules/local/purple/main'

include { CHANNEL_INPUTS_PURPLE } from './channel_inputs_purple'

workflow PURPLE_CALLING {
    take:
        // Sample data
        ch_inputs                    // channel: [mandatory] [ meta ]
        ch_amber                     // channel: [mandatory] [ meta, amber_dir ]
        ch_cobalt                    // channel: [mandatory] [ meta, cobalt_dir ]
        ch_smlv_somatic              // channel: [optional]  [ meta, pave_vcf ]
        ch_smlv_germline             // channel: [optional]  [ meta, pave_vcf ]
        ch_sv_somatic                // channel: [optional]  [ meta, gripss_vcf, gripss_tbi ]
        ch_sv_germline               // channel: [optional]  [ meta, gripss_vcf, gripss_tbi ]
        ch_sv_somatic_unfiltered     // channel: [optional]  [ meta, gripss_vcf, gripss_tbi ]

        // Reference data
        genome_fasta                 // channel: [mandatory] /path/to/genome_fasta
        genome_version               // channel: [mandatory] genome version
        genome_fai                   // channel: [mandatory] /path/to/genome_fai
        genome_dict                  // channel: [mandatory] /path/to/genome_dict
        gc_profile                   // channel: [mandatory] /path/to/gc_profile
        sage_known_hotspots_somatic  // channel: [mandatory] /path/to/sage_known_hotspots_somatic
        sage_known_hotspots_germline // channel: [mandatory] /path/to/sage_known_hotspots_germline
        driver_gene_panel            // channel: [mandatory] /path/to/driver_gene_panel
        ensembl_data_resources       // channel: [mandatory] /path/to/ensembl_data_resources/
        purple_germline_del          // channel: [mandatory] /path/to/purple_germline_del
        target_region_bed            // channel: [mandatory] /path/to/target_region_bed
        target_region_ratios         // channel: [mandatory] /path/to/target_region_ratios
        target_region_msi_indels     // channel: [mandatory] /path/to/target_region_msi_indels

        // Params
        run_config                   // channel: [mandatory] run configuration

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [ meta, amber_dir, cobalt_dir, sv_somatic_vcf, sv_somatic_tbi, sv_somatic_unfiltered_vcf, sv_somatic_unfiltered_tbi, sv_germline_vcf, sv_germline_tbi, smlv_somatic_vcf, smlv_germline_vcf ]
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
        // channel: [ meta_purple, amber_dir, cobalt_dir, sv_somatic_vcf, sv_somatic_tbi, sv_somatic_unfiltered_vcf, sv_somatic_unfiltered_tbi, sv_germline_vcf, sv_germline_tbi, smlv_somatic_vcf, smlv_germline_vcf ]
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
            genome_fasta,
            genome_version,
            genome_fai,
            genome_dict,
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
        purple_dir = ch_outputs  // channel: [ meta, purple_dir ]

        versions   = ch_versions // channel: [ versions.yml ]
}
