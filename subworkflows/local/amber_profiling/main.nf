//
// AMBER determines b-allele frequencies at predetermined positions
//

include { AMBER } from '../../../modules/local/amber/main'

workflow AMBER_PROFILING {
    take:
    // Sample data
    ch_inputs            // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor   // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal  // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_donor   // channel: [mandatory] [ meta, redux_dir ]

    // Reference data
    genome_fasta         // channel: [mandatory] /path/to/genome_fasta
    genome_version       // channel: [mandatory] genome version
    genome_fai           // channel: [mandatory] /path/to/genome_fai
    heterozygous_sites   // channel: [optional]  /path/to/heterozygous_sites
    target_regions_bed   // channel: [optional]  /path/to/target_regions_bed
    tumor_min_depth      // integer: [optional]  -tumor_min_depth argument value

    // Params
    sequencing_platform  // string:  [mandatory] sequencing platform
    purity_estimate_mode // boolean: [mandatory] Set purity estimate mode

    main:
    // Select input sources then sort
    // channel: runnable: [ meta, tumor_aln, tumor_idx, normal_aln, normal_idx ]
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

            def (tumor_aln, tumor_idx) = Utils.getTumorReduxDirAlignment(meta, redux_dir_tumor_selected)
            def (normal_aln, normal_idx) = Utils.getNormalReduxDirAlignment(meta, redux_dir_normal_selected)
            def (donor_aln, donor_idx) = Utils.getDonorReduxDirAlignment(meta, redux_dir_donor_selected)

            return [meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx]

        }
        .branch { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx ->

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.AMBER_DIR)
            def runnable_standard = ! purity_estimate_mode && tumor_aln && ! has_existing

            // TODO(SW): must improve handling through separation of sample information in meta; currently unable to provide ccfDNA AMBER directory in samplesheet
            def runnable_purity_estimate = purity_estimate_mode && normal_aln

            runnable: runnable_standard || runnable_purity_estimate
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_amber, tumor_aln, normal_aln, donor_aln, tumor_idx, normal_idx, donor_idx ]
    ch_amber_inputs = ch_inputs_sorted.runnable
        .map { meta, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx ->

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

            if (normal_aln) {
                meta_amber.normal_id = Utils.getNormalDnaSampleName(meta)
            }

            if (donor_aln) {
                meta_amber.donor_id = Utils.getDonorDnaSampleName(meta)
            }

            return [meta_amber, tumor_aln, tumor_idx, normal_aln, normal_idx, donor_aln, donor_idx]
        }

    // Run process
    AMBER(
        ch_amber_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        heterozygous_sites,
        target_regions_bed,
        tumor_min_depth,
        sequencing_platform,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, amber_dir ]
    ch_outputs = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('amber_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    amber_dir = ch_outputs // channel: [ meta, amber_dir ]
}
