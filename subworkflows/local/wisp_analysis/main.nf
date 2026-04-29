//
// WISP estimates tumor purity in longitudinal samples using WGS data of the primary
//

include { WISP } from '../../../modules/local/wisp/main'

workflow WISP_ANALYSIS {
    take:
    // Sample data
    ch_inputs                  // channel: [mandatory] [ meta ]
    ch_amber_out               // channel: [mandatory] [ meta, amber_dir ]
    ch_cobalt_out              // channel: [mandatory] [ meta, cobalt_dir ]
    ch_sage_somatic_append_out // channel: [mandatory] [ meta, sage_append_dir ]

    // Reference data
    genome_fasta               // channel: [mandatory] /path/to/genome_fasta
    genome_fai                 // channel: [mandatory] /path/to/genome_fai

    // Params
    targeted_mode              // boolean: [mandatory] Set targeted mode

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources and sort
    // channel: runnable: [ meta, ... ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = channels.WorkflowChannels.groupByMeta(
        ch_amber_out,
        ch_cobalt_out,
        ch_sage_somatic_append_out,
    )
        .branch { meta, longitudinal_amber_dir, longitudinal_cobalt_dir, longitudinal_sage_append_dir ->

            def primary_purple_dir = sample.Inputs.get(meta, sample.FileKey.PURPLE_DIR)
            def primary_amber_dir = sample.Inputs.get(meta, sample.FileKey.AMBER_DIR)
            def primary_normal_bam = sample.Inputs.get(meta, sample.FileKey.BAM_REDUX_DNA_NORMAL)

            def runnable
            if (targeted_mode) {
                runnable =
                    primary_purple_dir &&
                    longitudinal_sage_append_dir
            } else {
                runnable =
                    primary_purple_dir &&
                    longitudinal_sage_append_dir &&
                    longitudinal_cobalt_dir
            }

            def inputs = [:]
            inputs.meta = meta
            inputs.primary_purple_dir = primary_purple_dir
            inputs.primary_amber_dir = primary_amber_dir
            inputs.primary_normal_bam = primary_normal_bam
            inputs.longitudinal_amber_dir = longitudinal_amber_dir
            inputs.longitudinal_cobalt_dir = longitudinal_cobalt_dir
            inputs.longitudinal_sage_append_dir = longitudinal_sage_append_dir

            runnable: runnable
                return inputs
            skip: true
                return inputs.meta
        }

    // Create process input channel
    // channel: [ meta_wisp, ... ]
    ch_wisp_inputs = ch_inputs_sorted.runnable

        .map { inputs ->

            def meta = inputs.meta

            def meta_wisp = [
                key: meta.group_id,
                id: meta.group_id,
                subject_id: meta.subject_id,
                primary_id: sample.Inputs.getTumorDnaSampleNamePrimary(meta),
                longitudinal_id: sample.Inputs.getTumorDnaSampleNameLongitudinal(meta),
            ]

            inputs.meta = meta_wisp

            return inputs.values()
        }


    // Run process
    WISP(
        ch_wisp_inputs,
        genome_fasta,
        genome_fai,
        targeted_mode,
    )

    ch_versions = ch_versions.mix(WISP.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
