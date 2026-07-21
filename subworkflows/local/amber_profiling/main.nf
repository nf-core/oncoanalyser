//
// AMBER determines b-allele frequencies at predetermined positions
//

import Constants
import Utils

include { AMBER } from '../../../modules/local/amber/main'

workflow AMBER_PROFILING {
    take:
    // Sample data
    ch_inputs            // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor   // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal  // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_donor   // channel: [mandatory] [ meta, redux_dir ]

    // Reference data
    genome_version       // channel: [mandatory] genome version
    heterozygous_sites   // channel: [optional]  /path/to/heterozygous_sites
    target_regions_bed   // channel: [optional]  /path/to/target_regions_bed
    tumor_min_depth      // integer: [optional]  -tumor_min_depth argument value

    // Params
    sequencing_platform  // string:  [mandatory] sequencing platform
    purity_estimate_mode // boolean: [mandatory] Set purity estimate mode

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources then sort
    // channel: runnable: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_redux_dir_tumor,
        ch_redux_dir_normal,
        ch_redux_dir_donor,
    )
        .map { meta, redux_dir_tumor, redux_dir_normal, redux_dir_donor ->

            def redux_dir_tumor_selected = Utils.selectCurrentOrExisting(redux_dir_tumor, meta, Constants.INPUT.REDUX_DIR_TUMOR)
            def redux_dir_normal_selected = Utils.selectCurrentOrExisting(redux_dir_normal, meta, Constants.INPUT.REDUX_DIR_NORMAL)
            def redux_dir_donor_selected = Utils.selectCurrentOrExisting(redux_dir_donor, meta, Constants.INPUT.REDUX_DIR_DONOR)

            def (tumor_bam, tumor_bai) = Utils.getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)
            def (normal_bam, normal_bai) = Utils.getNormalReduxDirAlignment(meta, redux_dir_normal_selected)
            def (donor_bam, donor_bai) = Utils.getDonorReduxDirAlignment(meta, redux_dir_donor_selected)

            return [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai]

        }
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai ->

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.AMBER_DIR)
            def runnable_standard = ! purity_estimate_mode && tumor_bam && ! has_existing

            // TODO(SW): must improve handling through separation of sample information in meta; currently unable to provide ccfDNA AMBER directory in samplesheet
            def runnable_purity_estimate = purity_estimate_mode && normal_bam

            runnable: runnable_standard || runnable_purity_estimate
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_amber, tumor_bam, normal_bam, donor_bam, tumor_bai, normal_bai, donor_bai ]
    ch_amber_inputs = ch_inputs_sorted.runnable
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai ->

            def tumor_id
            if (purity_estimate_mode) {
                tumor_id = Utils.getTumorDnaSampleName(meta, primary: false)
            } else {
                tumor_id = Utils.getTumorDnaSampleName(meta, primary: true)
            }

            def meta_amber = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: tumor_id,
            ]

            if (normal_bam) {
                meta_amber.normal_id = Utils.getNormalDnaSampleName(meta)
            }

            if (donor_bam) {
                meta_amber.donor_id = Utils.getDonorDnaSampleName(meta)
            }

            return [meta_amber, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai]
        }

    // Run process
    AMBER(
        ch_amber_inputs,
        genome_version,
        heterozygous_sites,
        target_regions_bed,
        tumor_min_depth,
        sequencing_platform,
    )

    ch_versions = ch_versions.mix(AMBER.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, amber_dir ]
    ch_outputs = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(AMBER.out.amber_dir, ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    amber_dir = ch_outputs  // channel: [ meta, amber_dir ]

    versions  = ch_versions // channel: [ versions.yml ]
}
