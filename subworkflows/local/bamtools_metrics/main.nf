//
// Bam Tools calculates summary statistics for BAMs
//

include { BAMTOOLS } from '../../../modules/local/bamtools/main'

workflow BAMTOOLS_METRICS {
    take:
    // Sample data
    ch_inputs              // channel: [mandatory] [ meta ]
    ch_redux_dir_tumor     // channel: [mandatory] [ meta, redux_dir ]
    ch_redux_dir_normal    // channel: [mandatory] [ meta, redux_dir ]

    // Reference data
    genome_fasta           // channel: [mandatory] /path/to/genome_fasta
    genome_version         // channel: [mandatory] genome version
    genome_fai             // channel: [mandatory] /path/to/genome_fai
    driver_gene_panel      // channel: [mandatory] /path/to/driver_gene_panel
    ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/
    target_regions_bed     // channel: [optional]  /path/to/target_regions_bed

    main:
    // Select input sources then sort
    // channel: runnable: [ meta, aln, idx ]
    // channel: skip: [ meta ]
    ch_inputs_tumor_sorted = ch_redux_dir_tumor
        .map { meta, redux_dir ->

            def redux_dir_selected = Utils.selectCurrentOrExisting(redux_dir, meta, Constants.INPUT.REDUX_DIR_TUMOR)
            def (aln, idx) = Utils.getTumorReduxDirAlignment(meta, redux_dir_selected)

            return [meta, aln, idx]

        }
        .branch { meta, aln, idx ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAMTOOLS_DIR_TUMOR)
            runnable: aln && ! has_existing
            skip: true
                return meta
        }

    // channel: runnable: [ meta, aln, idx ]
    // channel: skip: [ meta ]
    ch_inputs_normal_sorted = ch_redux_dir_normal
        .map { meta, redux_dir ->

            def redux_dir_selected = Utils.selectCurrentOrExisting(redux_dir, meta, Constants.INPUT.REDUX_DIR_NORMAL)
            def (aln, idx) = Utils.getNormalReduxDirAlignment(meta, redux_dir_selected)

            return [meta, aln, idx]
        }
        .branch { meta, aln, idx ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAMTOOLS_DIR_NORMAL)
            runnable: aln && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_bamtools, aln, idx ]
    ch_bamtools_inputs = channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable.map { meta, aln, idx -> [meta, Utils.getTumorDnaSample(meta), 'tumor', aln, idx] },
            ch_inputs_normal_sorted.runnable.map { meta, aln, idx -> [meta, Utils.getNormalDnaSample(meta), 'normal', aln, idx] },
        )
        .map { meta, meta_sample, sample_type, aln, idx ->

            def meta_bamtools = [
                key: meta.group_id,
                id: "${meta.group_id}_${meta_sample.sample_id}",
                sample_id: meta_sample.sample_id,
                sample_type: sample_type,
            ]

            return [meta_bamtools, aln, idx]
        }

    // Run process
    BAMTOOLS(
        ch_bamtools_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        driver_gene_panel,
        ensembl_data_resources,
        target_regions_bed,
    )

    // Sort into a tumor and normal channel
    ch_bamtools_out = channel.topic('bamtools_metrics_dir')
        .branch { meta_bamtools, bamtools_dir ->
            assert ['tumor', 'normal'].contains(meta_bamtools.sample_type)
            tumor: meta_bamtools.sample_type == 'tumor'
            normal: meta_bamtools.sample_type == 'normal'
            placeholder: true
        }

    // Set outputs, restoring original meta
    // channel: [ meta, bamtools_dir ]
    ch_tumor_out = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_bamtools_out.tumor, ch_inputs),
            ch_inputs_tumor_sorted.skip.map { meta -> [meta, []] },
        )

    // channel: [ meta, bamtools_dir ]
    ch_normal_out = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_bamtools_out.normal, ch_inputs),
            ch_inputs_normal_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    tumor_dir  = ch_tumor_out  // channel: [ meta, bamtools_dir ]
    normal_dir = ch_normal_out // channel: [ meta, bamtools_dir ]
}
