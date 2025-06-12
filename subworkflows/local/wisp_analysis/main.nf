//
// WISP estimates tumor purity in longitudinal samples using WGS data of the primary
//

import Constants
import Utils

include { WISP } from '../../../modules/local/wisp/main'

workflow WISP_ANALYSIS {
    take:
    // Sample data
    ch_inputs                  // channel: [mandatory] [ meta ]
    ch_amber_out               // channel: [mandatory] [ meta, amber_dir ]
    ch_cobalt_out              // channel: [mandatory] [ meta, cobalt_dir ]
    ch_sage_somatic_append_out // channel: [mandatory] [ meta, sage_append_dir ]

    // Reference data
    genome_fasta     // channel: [mandatory] /path/to/genome_fasta
    genome_fai       // channel: [mandatory] /path/to/genome_fai

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    ch_primary_amber_dir = ch_inputs
        .map { meta ->
            return [meta, Utils.selectCurrentOrExisting([], meta, Constants.INPUT.AMBER_DIR)]
        }

    ch_primary_purple_dir = ch_inputs
        .map { meta ->
            return [meta, Utils.selectCurrentOrExisting([], meta, Constants.INPUT.PURPLE_DIR)]
        }

    // Select input sources and sort
    // channel: runnable: [ meta, ... ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_sage_somatic_append_out,
        ch_amber_out,
        ch_cobalt_out,
        ch_primary_amber_dir,
        ch_primary_purple_dir,
    )
        .branch { meta, sage_append_dir, amber_dir, cobalt_dir, primary_amber_dir, primary_purple_dir ->

            def purity_estimate_mode = Utils.getEnumFromString(params.purity_estimate_mode, Constants.RunMode)

            def runnable = purity_estimate_mode === Constants.RunMode.WGTS
                ? primary_purple_dir && primary_amber_dir && sage_append_dir && amber_dir && cobalt_dir
                : primary_purple_dir                      && sage_append_dir

            runnable: runnable
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_wisp, ... ]
    ch_wisp_inputs = ch_inputs_sorted.runnable

        .map { meta, sage_append_dir, amber_dir, cobalt_dir, primary_amber_dir, primary_purple_dir ->

            def tumor_dna_id = Utils.getTumorDnaSampleName(meta, primary: true)

            def meta_wisp = [
                key: meta.group_id,
                id: meta.group_id,
                subject_id: meta.subject_id,
                primary_id: Utils.getTumorDnaSampleName(meta, primary: true),
                sample_id: Utils.getTumorDnaSampleName(meta),
            ]

            return [meta_wisp, sage_append_dir, amber_dir, cobalt_dir, primary_amber_dir, primary_purple_dir]
        }


    // Run process
    WISP(
        ch_wisp_inputs,
        genome_fasta,
        genome_fai,
    )

    ch_versions = ch_versions.mix(WISP.out.versions)

    emit:
    versions     = ch_versions     // channel: [ versions.yml ]
}
