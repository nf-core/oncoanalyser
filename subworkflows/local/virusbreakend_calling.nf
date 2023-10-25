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

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Sort inputs
        // NOTE(SW): VIRUSBreakend inputs are not allowed in the samplesheet, so aren't considered
        // channel: [ meta ]
        ch_inputs_sorted = ch_inputs
            .branch { meta ->

                def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.VIRUSINTERPRETER_DIR)

                runnable: Utils.hasTumorDnaBam(meta) && !has_existing
                skip: true
            }

        //
        // MODULE: VIRUSBreakend
        //
        // Create process input channel
        // channel: [ meta_virus, tumor_bam ]
        ch_virusbreakend_inputs = ch_inputs_sorted.runnable
            .map { meta ->

                def meta_virus = [
                    key: meta.group_id,
                    id: meta.group_id,
                    sample_id: Utils.getTumorDnaSampleName(meta),
                ]

                return [meta_virus, Utils.getTumorDnaBam(meta)]
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

        ch_versions = ch_versions.mix(VIRUSBREAKEND.out.versions)

        //
        // MODULE: Virus Interpreter
        //
        // Select input sources
        // channel: [ meta, virus_tsv, purple_dir, metrics ]
        ch_virusinterpreter_inputs_selected = WorkflowOncoanalyser.groupByMeta(
              WorkflowOncoanalyser.restoreMeta(VIRUSBREAKEND.out.tsv, ch_inputs),
              ch_purple,
              ch_bamtools_somatic,
        )
            .map { meta, virus_tsv, purple_dir, metrics ->

                def inputs = [
                    virus_tsv,
                    Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
                    Utils.selectCurrentOrExisting(metrics, meta, Constants.INPUT.BAMTOOLS_TUMOR),
                ]

                return [meta, *inputs]
            }

        // Sort inputs
        // channel: [ meta, virus_tsv, purple_dir, metrics ]
        // channel: skip: [ meta ]
        ch_virusinterpreter_inputs_sorted = ch_virusinterpreter_inputs_selected
            .branch { meta, virus_tsv, purple_dir, metrics ->
                runnable: virus_tsv && purple_dir && metrics
                skip: true
                    return meta
            }

        // Create process input channel
        // channel: [ meta_virus, virus_tsv, purple_dir, metrics ]
        ch_virusinterpreter_inputs = ch_virusinterpreter_inputs_sorted.runnable
            .map { d ->

                def meta = d[0]
                def inputs = d[1..-1]

                def meta_virus = [
                    key: meta.group_id,
                    id: meta.group_id,
                    sample_id: Utils.getTumorDnaSampleName(meta),
                ]

                return [meta_virus, *inputs]
            }

        // Run process
        VIRUSINTERPRETER(
            ch_virusinterpreter_inputs,
            virus_taxonomy_db,
            virus_reporting_db,
        )

        ch_versions = ch_versions.mix(VIRUSINTERPRETER.out.versions)

        // Set outputs, restoring original meta
        // channel: [ meta, virusinterpreter_dir ]
        ch_outputs = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(VIRUSINTERPRETER.out.virusinterpreter_dir, ch_inputs),
                ch_virusinterpreter_inputs_sorted.skip.map { meta -> [meta, []] },
                ch_inputs_sorted.skip.map { meta -> [meta, []] },
            )

    emit:
        virusinterpreter_dir = ch_outputs  // channel: [ meta, virusinterpreter_dir ]

        versions             = ch_versions // channel: [ versions.yml ]
}
