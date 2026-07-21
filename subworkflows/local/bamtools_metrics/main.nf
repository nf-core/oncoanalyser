//
// Bam Tools calculates summary statistics for BAMs
//

import Constants
import Utils

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
    driver_gene_panel      // channel: [mandatory] /path/to/driver_gene_panel
    ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/
    target_regions_bed     // channel: [optional]  /path/to/target_regions_bed

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources then sort
    // channel: runnable: [ meta, bam, bai ]
    // channel: skip: [ meta ]
    ch_inputs_tumor_sorted = ch_redux_dir_tumor
        .map { meta, redux_dir ->

            def redux_dir_selected = Utils.selectCurrentOrExisting(redux_dir, meta, Constants.INPUT.REDUX_DIR_TUMOR)
            def (bam, bai) = Utils.getTumorReduxDirAlignment(meta, redux_dir_selected)

            return [meta, bam, bai]

        }
        .branch { meta, bam, bai ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAMTOOLS_DIR_TUMOR)
            runnable: bam && ! has_existing
            skip: true
                return meta
        }

    // channel: runnable: [ meta, bam, bai ]
    // channel: skip: [ meta ]
    ch_inputs_normal_sorted = ch_redux_dir_normal
        .map { meta, redux_dir ->

            def redux_dir_selected = Utils.selectCurrentOrExisting(redux_dir, meta, Constants.INPUT.REDUX_DIR_NORMAL)
            def (bam, bai) = Utils.getNormalReduxDirAlignment(meta, redux_dir_selected)

            return [meta, bam, bai]
        }
        .branch { meta, bam, bai ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAMTOOLS_DIR_NORMAL)
            runnable: bam && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_bamtools, bam, bai ]
    ch_bamtools_inputs = Channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable.map { meta, bam, bai -> [meta, Utils.getTumorDnaSample(meta), 'tumor', bam, bai] },
            ch_inputs_normal_sorted.runnable.map { meta, bam, bai -> [meta, Utils.getNormalDnaSample(meta), 'normal', bam, bai] },
        )
        .map { meta, meta_sample, sample_type, bam, bai ->

            def meta_bamtools = [
                key: meta.group_id,
                id: "${meta.group_id}_${meta_sample.sample_id}",
                sample_id: meta_sample.sample_id,
                sample_type: sample_type,
            ]

            return [meta_bamtools, bam, bai]
        }

    // Run process
    BAMTOOLS(
        ch_bamtools_inputs,
        genome_fasta,
        genome_version,
        driver_gene_panel,
        ensembl_data_resources,
        target_regions_bed,
    )

    ch_versions = ch_versions.mix(BAMTOOLS.out.versions)

    // Sort into a tumor and normal channel
    ch_bamtools_out = BAMTOOLS.out.bamtools_dir
        .branch { meta_bamtools, bamtools_dir ->
            assert ['tumor', 'normal'].contains(meta_bamtools.sample_type)
            tumor: meta_bamtools.sample_type == 'tumor'
            normal: meta_bamtools.sample_type == 'normal'
            placeholder: true
        }

    // Set outputs, restoring original meta
    // channel: [ meta, bamtools_dir ]
    ch_tumor_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_bamtools_out.tumor, ch_inputs),
            ch_inputs_tumor_sorted.skip.map { meta -> [meta, []] },
        )

    // channel: [ meta, bamtools_dir ]
    ch_normal_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_bamtools_out.normal, ch_inputs),
            ch_inputs_normal_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    tumor_dir  = ch_tumor_out  // channel: [ meta, bamtools_dir ]
    normal_dir = ch_normal_out // channel: [ meta, bamtools_dir ]

    versions   = ch_versions   // channel: [ versions.yml ]
}
