//
// SAGE is a precise and highly sensitive somatic SNV, MNV and small INDEL caller
//

include { SAGE_VISUALISER } from '../../../modules/local/sage/visualiser/main'

workflow SAGE_PLOTTING {
    take:
    // Sample data
    ch_inputs                    // channel: [mandatory] [ meta ]
    ch_tumor_bam                 // channel: [mandatory] [ meta, bam, bai ]
    ch_normal_bam                // channel: [mandatory] [ meta, bam, bai ]
    ch_donor_bam                 // channel: [mandatory] [ meta, bam, bai ]
    ch_tumor_dir                 // channel: [mandatory] [ meta, redux_dir ]
    ch_normal_dir                // channel: [mandatory] [ meta, redux_dir ]
    ch_donor_dir                 // channel: [mandatory] [ meta, redux_dir ]
    ch_purple                    // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    genome_fasta                 // channel: [mandatory] /path/to/genome_fasta
    genome_version               // channel: [mandatory] genome version
    genome_fai                   // channel: [mandatory] /path/to/genome_fai
    genome_dict                  // channel: [mandatory] /path/to/genome_dict
    sage_pon                     // channel: [mandatory] /path/to/sage_pon
    sage_known_hotspots_somatic  // channel: [mandatory] /path/to/sage_known_hotspots_somatic
    sage_highconf_regions        // channel: [mandatory] /path/to/sage_highconf_regions
    ensembl_data_resources       // channel: [mandatory] /path/to/ensembl_data_resources/

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Select input sources. Route inputs
    // channel: runnable: { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, [redux_tsv, ...], purple_dir }
    // channel: skip: [ meta ]
    ch_inputs_sorted = channels.WorkflowChannels.groupByMeta(
        flatten_mode: 'singletons_only',
        ch_tumor_bam, ch_normal_bam, ch_donor_bam,
        ch_tumor_dir, ch_normal_dir, ch_donor_dir,
        ch_purple,
    )
        .map { meta, tumor_bam_bai, normal_bam_bai, donor_bam_bai, tumor_dir, normal_dir, donor_dir, purple_dir ->

            def (tumor_bam, tumor_bai) = Inputs.resolveReduxBamBai(tumor_bam_bai, meta, samplesheet.SampleType.TUMOR)
            def (normal_bam, normal_bai) = Inputs.resolveReduxBamBai(normal_bam_bai, meta, samplesheet.SampleType.NORMAL)
            def (donor_bam, donor_bai) = Inputs.resolveReduxBamBai(donor_bam_bai, meta, samplesheet.SampleType.DONOR)

            def tumor_tsvs = Inputs.resolveReduxTsvFiles(tumor_dir, meta, samplesheet.SampleType.TUMOR)
            def normal_tsvs = Inputs.resolveReduxTsvFiles(normal_dir, meta, samplesheet.SampleType.NORMAL)
            def donor_tsvs = Inputs.resolveReduxTsvFiles(donor_dir, meta, samplesheet.SampleType.DONOR)
            def redux_tsvs = [ *tumor_tsvs, *normal_tsvs, *donor_tsvs ]

            purple_dir = Inputs.preferUserProvidedInput(purple_dir, meta, Inputs.KEY.PURPLE_DIR)

            def inputs = [
                meta: meta,
                tumor_bam: tumor_bam,
                tumor_bai: tumor_bai,
                normal_bam: normal_bam,
                normal_bai: normal_bai,
                donor_bam: donor_bam,
                donor_bai: donor_bai,
                redux_tsvs: redux_tsvs,
                purple_dir: purple_dir,
            ]

            return inputs
        }
        .branch { inputs ->
            runnable: inputs.tumor_bam && inputs.purple_dir
                return inputs
            skip: true
                return inputs.meta
        }

    // Create process input channel
    // channel: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, [redux_tsv, ...], purple_smlv_vcf ]
    ch_sage_plotting_inputs = ch_inputs_sorted.runnable
        .map { inputs ->

            def meta = inputs.meta

            def meta_sage = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Inputs.getTumorDnaSampleName(meta),
                normal_id: inputs.normal_bam ? Inputs.getNormalDnaSampleName(meta) : null,
                donor_id: inputs.donor_bam ? Inputs.getDonorDnaSampleName(meta) : null,
            ]

            def purple_smlv_vcf = Inputs.resolvePurpleSomaticVcf(inputs.purple_dir, meta)
            def purple_smlv_vcf_tbi = Inputs.resolvePurpleSomaticVcfTbi(inputs.purple_dir, meta)

            return [
                meta_sage,
                inputs.tumor_bam,
                inputs.normal_bam,
                inputs.donor_bam,
                inputs.tumor_bai,
                inputs.normal_bai,
                inputs.donor_bai,
                inputs.redux_tsvs,
                purple_smlv_vcf,
                purple_smlv_vcf_tbi
            ]
        }

    // Run process
    SAGE_VISUALISER(
        ch_sage_plotting_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        sage_pon,
        sage_known_hotspots_somatic,
        sage_highconf_regions,
        ensembl_data_resources,
    )

    ch_versions = ch_versions.mix(SAGE_VISUALISER.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, sage_dir ]
    ch_visualiser_dir_out = Channel.empty()
        .mix(
            channels.WorkflowChannels.restoreMeta(SAGE_VISUALISER.out.sage_vis_dir, ch_inputs),
            channels.PlaceholderChannels.toolDir(ch_inputs_sorted.skip),
        )

    emit:
    visualiser_dir  = ch_visualiser_dir_out // channel: [ meta, sage_dir ]
    versions        = ch_versions           // channel: [ versions.yml ]
}
