//
// MultiQC aggregates and collates metrics for QC review
//

include { MULTIQC } from '../../../modules/local/multiqc/main'

include { paramsSummaryMap } from 'plugin/nf-schema'

include { paramsSummaryMultiqc   } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../../subworkflows/local/utils_nfcore_oncoanalyser_pipeline'

workflow MULTIQC_REPORTING {
    take:
    // Sample data
    ch_bamtools_dir_tumor     // channel: [mandatory] [ meta, bamtools_dir ]
    ch_bamtools_dir_normal    // channel: [optional]  [ meta, bamtools_dir ]
    ch_amber_dir              // channel: [mandatory] [ meta, amber_dir ]
    ch_purple_dir             // channel: [mandatory] [ meta, purple_dir ]
    ch_align_rna_qc_tumor_out // channel: [mandatory] [ meta, star_log, rna_md_metrics ]

    // Other
    ch_collated_versions      // channel: [mandatory] [ collated_versions.yml ]
    custom_config             //  string: [optional]  Custom configuration for MultiQC
    custom_desc               //  string: [optional]  Custom methods description for MultiQC
    custom_logo               //  string: [optional]  Custom logo for MultiQC

    main:
    // Select input sources then sort
    // channel: [ meta, bamtools_tumor_dir, bamtools_normal_dir, amber_dir, purple_dir, star_log, rna_md_metrics ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_bamtools_dir_tumor,
        ch_bamtools_dir_normal,
        ch_amber_dir,
        ch_purple_dir,
        ch_align_rna_qc_tumor_out,
    )
        .map { meta, bamtools_dir_tumor, bamtools_dir_normal, amber_dir, purple_dir, star_log, rna_md_metrics ->

            // NOTE(SW): will not implement ability for user to provide RNA alignment QC metrics

            return [
                meta,
                Utils.selectCurrentOrExisting(bamtools_dir_tumor, meta, Constants.INPUT.BAMTOOLS_DIR_TUMOR),
                Utils.selectCurrentOrExisting(bamtools_dir_normal, meta, Constants.INPUT.BAMTOOLS_DIR_NORMAL),
                Utils.selectCurrentOrExisting(amber_dir, meta, Constants.INPUT.AMBER_DIR),
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
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
    // channel: [ group_ids, file_group_ids, files ]
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

        // channel: [ meta_multiqc, file_group_ids, files ]
        .map { metas_and_files  ->

            def files = []
            def files_group_ids = []
            def group_sample_ids = [:]

            metas_and_files.each { meta, files_nested ->

                files_nested.each { f ->
                    if (! f) {
                        return
                    }

                    files_group_ids << meta.group_id
                    files << f
                }

              // NOTE(SW): not handled here are cases with no input files, unlikely scenario (impossible?)

              group_sample_ids[meta.group_id] = [
                  'normal_dna_id': Utils.getNormalDnaSampleName(meta),
                  'tumor_dna_id': Utils.getTumorDnaSampleName(meta),
                  'tumor_rna_id': Utils.getTumorRnaSampleName(meta),
              ]

            }

            // Perform stable sort on aggregated data and without channels ops to ensure no NF language assumptions are broken
            // This sort is done to remove some non-deterministic behaviour in MultiQC resulting from unstable input ordering
            // NOTE(SW): nonetheless approaches to create other inputs are non-deterministic and hence MultiQC is never cached
            def (files_group_ids_sorted, files_sorted) = [files_group_ids, files]
                // Pair: [ [gid_e, f], [gid_a, f], ...]
                .transpose()
                // Sort: [ [gid_a, f], [gid_b, f], ...]
                // Sorting done on primarily on group_id (index 0) then secondarily on filename (index 1)
                // NOTE(SW): assumes common prefix in filenames from each output stage that is consistent between samples
                .sort { d -> d[1].name }
                .sort { d -> d[0] }
                // Separate: [ [gid_a, gid_b, ...], [f, f, ..] ]
                .transpose()

            return [group_sample_ids, files_group_ids_sorted, files_sorted]

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
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
    )

    // Set outputs
    ch_outputs = channel.topic('multiqc_report').toList()

    emit:
    report = ch_outputs // path: multiqc_report
}
