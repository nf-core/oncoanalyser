//
// VIRUSBreakend and Virus Interpreter identify viral content and insertion sites
//

include { VIRUSBREAKEND    } from '../../../modules/local/virusbreakend/main'
include { VIRUSINTERPRETER } from '../../../modules/local/virusinterpreter/main'

workflow VIRUSBREAKEND_CALLING {
    take:
    // Sample data
    ch_inputs           // channel: [mandatory] [ meta ]
    ch_tumor_bam        // channel: [mandatory] [ meta, bam, bai ]
    ch_purple           // channel: [mandatory] [ meta, purple_dir ]
    ch_bamtools_somatic // channel: [mandatory] [ meta, metrics ]

    // Reference data
    genome_fasta        // channel: [mandatory] /path/to/genome_fasta
    genome_fai          // channel: [mandatory] /path/to/genome_fai
    genome_dict         // channel: [mandatory] /path/to/genome_dict
    genome_gridss_index // channel: [mandatory] /path/to/genome_gridss_index
    virusbreakenddb     // channel: [mandatory] /path/to/virusbreakenddb/
    virus_taxonomy_db   // channel: [mandatory] /path/to/virus_taxonomy_db
    virus_reporting_db  // channel: [mandatory] /path/to/virus_reporting_db
    virus_blocklist_db  // channel: [mandatory] /path/to/virus_blocklist_db

    // Params
    gridss_config          // channel: [optional] /path/to/gridss_config

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Sort inputs
    // NOTE(SW): VIRUSBreakend inputs are not allowed in the samplesheet, so aren't considered
    // channel: [ meta, tumor_bam, tumor_bai ]
    ch_inputs_sorted = ch_tumor_bam
        .map { meta, tumor_bam, tumor_bai ->
            return [
                meta,
                Inputs.preferUserProvidedInput(tumor_bam, meta, SampleSheetFields.INPUT.BAM_REDUX_DNA_TUMOR),
                Inputs.preferUserProvidedInput(tumor_bai, meta, SampleSheetFields.INPUT.BAI_DNA_TUMOR),
            ]
        }
        .branch { meta, tumor_bam, tumor_bai ->
            def has_existing = Inputs.hasExistingInput(meta, SampleSheetFields.INPUT.VIRUSINTERPRETER_DIR)
            runnable: tumor_bam && !has_existing
            skip: true
                return meta
        }

    //
    // MODULE: VIRUSBreakend
    //
    // Create process input channel
    // channel: [ meta_virus, tumor_bam ]
    ch_virusbreakend_inputs = ch_inputs_sorted.runnable
        .map { meta, tumor_bam, tumor_bai ->

            def meta_virus = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Inputs.getTumorDnaSampleName(meta),
            ]

            return [meta_virus, tumor_bam]
        }

    // Run process
    VIRUSBREAKEND(
        ch_virusbreakend_inputs,
        genome_fasta,
        genome_fai,
        genome_dict,
        genome_gridss_index,
        virusbreakenddb,
        gridss_config,
    )

    ch_versions = ch_versions.mix(VIRUSBREAKEND.out.versions)

    //
    // MODULE: Virus Interpreter
    //
    // Select input sources
    // channel: [ meta, virus_tsv, purple_dir, metrics ]
    ch_virusinterpreter_inputs_selected = WorkflowChannels.groupByMeta(
        WorkflowChannels.restoreMeta(VIRUSBREAKEND.out.tsv, ch_inputs),
        ch_purple,
        ch_bamtools_somatic,
    )
        .map { meta, virus_tsv, purple_dir, somatic_metrics ->

            def inputs = [
                virus_tsv,
                Inputs.preferUserProvidedInput(purple_dir, meta, SampleSheetFields.INPUT.PURPLE_DIR),
                Inputs.preferUserProvidedInput(somatic_metrics, meta, SampleSheetFields.INPUT.BAMTOOLS_DIR_TUMOR),
            ]

            return [meta, *inputs]
        }

    // Sort inputs
    // channel: [ meta, virus_tsv, purple_dir, metrics ]
    // channel: skip: [ meta ]
    ch_virusinterpreter_inputs_sorted = ch_virusinterpreter_inputs_selected
        .branch { meta, virus_tsv, purple_dir, somatic_metrics ->
            runnable: virus_tsv && purple_dir && somatic_metrics
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
                sample_id: Inputs.getTumorDnaSampleName(meta),
            ]

            return [meta_virus, *inputs]
        }

    // Run process
    VIRUSINTERPRETER(
        ch_virusinterpreter_inputs,
        virus_taxonomy_db,
        virus_reporting_db,
        virus_blocklist_db,
    )

    ch_versions = ch_versions.mix(VIRUSINTERPRETER.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, virusinterpreter_dir ]
    ch_outputs = Channel.empty()
        .mix(
            WorkflowChannels.restoreMeta(VIRUSINTERPRETER.out.virusinterpreter_dir, ch_inputs),
            PlaceholderChannels.toolDir(ch_virusinterpreter_inputs_sorted.skip),
            PlaceholderChannels.toolDir(ch_inputs_sorted.skip),
        )

    emit:
    virusinterpreter_dir = ch_outputs  // channel: [ meta, virusinterpreter_dir ]

    versions             = ch_versions // channel: [ versions.yml ]
}
