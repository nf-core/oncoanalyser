//
// MultiQC aggregates and collates metrics for QC review
//

nextflow.enable.types = true

include { MULTIQC  } from '../../../modules/local/multiqc/main'

include { paramsSummaryMap  } from 'plugin/nf-schema'

include { paramsSummaryMultiqc  } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { FileType                 } from '../utils_nfcore_oncoanalyser_pipeline/types_enums'
include { methodsDescriptionText  } from '../../../subworkflows/local/utils_nfcore_oncoanalyser_pipeline'
include { groupByMeta              } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { joinMeta                 } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { restoreMeta              } from '../utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { getInput                 } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalDnaSample       } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getNormalDnaSampleName   } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSample        } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorDnaSampleName    } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { getTumorRnaSampleName    } from '../utils_nfcore_oncoanalyser_pipeline/accessors_samples'
include { selectCurrentOrExisting  } from '../utils_nfcore_oncoanalyser_pipeline/utils'

workflow MULTIQC_REPORTING {
    take:
    // Sample data
    ch_bamtools_dir_tumor: Channel<Tuple<Map, Path>>     // channel: [mandatory] [ meta, bamtools_dir ]
    ch_bamtools_dir_normal: Channel<Tuple<Map, Path>>    // channel: [optional]  [ meta, bamtools_dir ]
    ch_amber_dir: Channel<Tuple<Map, Path>>              // channel: [mandatory] [ meta, amber_dir ]
    ch_purple_dir: Channel<Tuple<Map, Path>>             // channel: [mandatory] [ meta, purple_dir ]
    ch_align_rna_qc_tumor_out: Channel<Tuple<Map, Path, Path>> // channel: [mandatory] [ meta, star_log, rna_md_metrics ]

    // Other
    ch_collated_versions: Channel<Path>      // channel: [mandatory] [ collated_versions.yml ]
    custom_config: String?             //  string: [optional]  Custom configuration for MultiQC
    custom_desc: String?               //  string: [optional]  Custom methods description for MultiQC
    custom_logo: String?               //  string: [optional]  Custom logo for MultiQC

    main:
    // Select input sources then sort
    // channel: [ meta, bamtools_tumor_dir, bamtools_normal_dir, amber_dir, purple_dir, star_log, rna_md_metrics ]
    ch_inputs_sorted = groupByMeta([
        ch_bamtools_dir_tumor,
        ch_bamtools_dir_normal,
        ch_amber_dir,
        ch_purple_dir,
        ch_align_rna_qc_tumor_out,
    ])
        .map { meta, bamtools_dir_tumor, bamtools_dir_normal, amber_dir, purple_dir, star_log, rna_md_metrics ->

            // NOTE(SW): will not implement ability for user to provide RNA alignment QC metrics

            return [
                meta,
                selectCurrentOrExisting(bamtools_dir_tumor, getInput(getTumorDnaSample(meta), FileType.BAMTOOLS_DIR)),
                selectCurrentOrExisting(bamtools_dir_normal, getInput(getNormalDnaSample(meta), FileType.BAMTOOLS_DIR)),
                selectCurrentOrExisting(amber_dir, getInput(getTumorDnaSample(meta), FileType.AMBER_DIR)),
                selectCurrentOrExisting(purple_dir, getInput(getTumorDnaSample(meta), FileType.PURPLE_DIR)),
                star_log,
                rna_md_metrics,
            ]

        }
        .branch { meta, bamtools_dir_tumor, bamtools_dir_normal, amber_dir, purple_dir, star_log, rna_md_metrics ->

            runnable: bamtools_dir_tumor || bamtools_dir_normal || amber_dir || purple_dir || star_log || rna_md_metrics
            skip: true
                return meta
        }

    // Create linkage between identifiers and output files
    // channel: [ case_ids, file_case_ids, files ]
    ch_multiqc_files_sample = ch_inputs_sorted.runnable

        // channel: [ meta, [ file_one, file_two, ... ] ]
        .map { d ->
            def meta = d[0]
            def fs = d[1..-1]
            return [meta, fs]
        }

        // channel: [
        //     [
        //         [ meta, [ file_one, file_two, ... ],  // group_a
        //         [ meta, [ file_one, file_two, ... ],  // group_b
        //         ...
        //     ]
        // ]
        .collect(flat: false)

        // channel: [ meta_multiqc, file_case_ids, files ]
        .map { metas_and_files  ->

            def files = []
            def files_case_ids = []
            def group_sample_ids = [:]

            metas_and_files.each { meta, files_nested ->

                files_nested.each { f ->
                    if (! f) {
                        return
                    }

                    files_case_ids << meta.case_id
                    files << f
                }

              // NOTE(SW): not handled here are cases with no input files, unlikely scenario (impossible?)

              group_sample_ids[meta.case_id] = [
                  'normal_dna_id': getNormalDnaSampleName(meta),
                  'tumor_dna_id': getTumorDnaSampleName(meta),
                  'tumor_rna_id': getTumorRnaSampleName(meta),
              ]

            }

            // Perform stable sort on aggregated data and without channels ops to ensure no NF language assumptions are broken
            // This sort is done to remove some non-deterministic behaviour in MultiQC resulting from unstable input ordering
            // NOTE(SW): nonetheless approaches to create other inputs are non-deterministic and hence MultiQC is never cached
            def (files_case_ids_sorted, files_sorted) = [files_case_ids, files]
                // Pair: [ [gid_e, f], [gid_a, f], ...]
                .transpose()
                // Sort: [ [gid_a, f], [gid_b, f], ...]
                // Sorting done on primarily on case_id (index 0) then secondarily on filename (index 1)
                // NOTE(SW): assumes common prefix in filenames from each output stage that is consistent between samples
                .sort { d -> d[1].name }
                .sort { d -> d[0] }
                // Separate: [ [gid_a, gid_b, ...], [f, f, ..] ]
                .transpose()

            return [group_sample_ids, files_case_ids_sorted, files_sorted]

        }

    // Channel for other input files; essentially unused but retained for nf-core
    ch_multiqc_files = channel.empty()

    // nf-core boilerplate
    ch_multiqc_config = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = custom_config ? channel.fromPath(custom_config, checkIfExists: true) : channel.empty()
    ch_multiqc_logo = custom_logo ? channel.fromPath(custom_logo, checkIfExists: true) : channel.empty()

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = custom_desc ? file(custom_desc, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    // Run process
    MULTIQC(
        ch_multiqc_files_sample,
        ch_multiqc_files.toList(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
    )

    // Set outputs
    ch_outputs = channel.topic('multiqc_report').toList()

    emit:
    report = ch_outputs // path: multiqc_report
}
