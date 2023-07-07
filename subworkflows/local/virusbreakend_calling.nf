//
// VIRUSBreakend and Virus Interpreter identify viral content and insertion sites
//

import Constants
import Utils

include { VIRUSBREAKEND    } from '../../modules/local/virusbreakend/main'
include { VIRUSINTERPRETER } from '../../modules/local/virusinterpreter/main'

workflow VIRUSBREAKEND_CALLING {
    take:
        // Sample data
        ch_inputs              // channel: [mandatory] [ meta ]
        ch_purple              // channel: [mandatory] [ meta, purple_dir ]
        ch_bamtools_somatic    // channel: [mandatory] [ meta, metrics ]

        // Reference data
        genome_fasta           // channel: [mandatory] /path/to/genome_fasta
        genome_fai             // channel: [mandatory] /path/to/genome_fai
        genome_dict            // channel: [mandatory] /path/to/genome_dict
        genome_bwa_index       // channel: [mandatory] /path/to/genome_bwa_index/
        genome_bwa_index_image // channel: [mandatory] /path/to/genome_bwa_index_image
        genome_gridss_index    // channel: [mandatory] /path/to/genome_gridss_index
        virusbreakenddb        // channel: [mandatory] /path/to/virusbreakenddb/
        virus_taxonomy_db      // channel: [mandatory] /path/to/virus_taxonomy_db
        virus_reporting_db     // channel: [mandatory] /path/to/virus_reporting_db

        // Params
        gridss_config          // channel: [optional] /path/to/gridss_config
        run_config             // channel: [mandatory] run configuration

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // VIRUSBreakend
        // Create inputs and create process-specific meta
        // channel: [ meta_virus, tumor_bam ]
        ch_virusbreakend_inputs = ch_inputs
            .map { meta ->
                def meta_virus = [
                    key: meta.id,
                    id: meta.id,
                ]
                return [meta_virus, Utils.getTumorBam(meta, run_config.mode)]
            }

        // Run process
        VIRUSBREAKEND(
            ch_virusbreakend_inputs,
            gridss_config,
            genome_fasta,
            genome_fai,
            genome_dict,
            genome_bwa_index,
            genome_bwa_index_image,
            genome_gridss_index,
            virusbreakenddb,
        )

        // Create inputs and create process-specific meta
        if (run_config.stages.purple) {
            ch_virusinterpreter_inputs_purple = ch_purple
        } else {
            ch_virusinterpreter_inputs_purple = WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PURPLE_DIR)
        }

        // channel: [ meta, purple_qc, wgs_metrics ]
        ch_virusinterpreter_inputs_purple_files = ch_virusinterpreter_inputs_purple
            .filter { it[0] != Constants.PLACEHOLDER_META }
            .map { meta, purple_dir ->
                def tumor_id = Utils.getTumorSampleName(meta, run_config.mode)
                def purple_purity = file(purple_dir).resolve("${tumor_id}.purple.purity.tsv")
                def purple_qc = file(purple_dir).resolve("${tumor_id}.purple.qc")

                // Require both purity and QC files from the PURPLE directory
                if (!purple_purity.exists() || !purple_qc.exists()) {
                    return Constants.PLACEHOLDER_META
                }
                return [meta, purple_purity, purple_qc]
            }
            .filter { it != Constants.PLACEHOLDER_META }

        // channel: [ meta, virus_tsv, purple_purity, purple_qc, wgs_metrics ]
        ch_virusinterpreter_inputs_full = WorkflowOncoanalyser.groupByMeta(
            WorkflowOncoanalyser.restoreMeta(VIRUSBREAKEND.out.tsv, ch_inputs),
            ch_virusinterpreter_inputs_purple_files,
            run_config.stages.bamtools ? ch_bamtools_somatic : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.INPUT_BAMTOOLS_TUMOR),
        )

        // Virus Interpreter
        // Create inputs and create process-specific meta
        // channel: [ meta_virus, virus_tsv, purple_purity, purple_qc, wgs_metrics ]
        ch_virusinterpreter_inputs = ch_virusinterpreter_inputs_full
            .map {
                def meta = it[0]
                def meta_virus = [
                    key: meta.id,
                    id: Utils.getTumorSampleName(meta, run_config.mode),
                ]
                return [meta_virus, *it[1..-1]]
            }

        // Run process
        VIRUSINTERPRETER(
            ch_virusinterpreter_inputs,
            virus_taxonomy_db,
            virus_reporting_db,
        )

        // Set outputs, restoring original meta
        ch_outputs = WorkflowOncoanalyser.restoreMeta(VIRUSINTERPRETER.out.virusinterpreter, ch_inputs)
        ch_versions = ch_versions.mix(
            VIRUSINTERPRETER.out.versions,
            VIRUSBREAKEND.out.versions,
        )

    emit:
        virusinterpreter = ch_outputs  // channel: [ meta, virusinterpreter ]

        versions         = ch_versions // channel: [ versions.yml ]
}
