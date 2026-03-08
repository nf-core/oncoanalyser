//
// SAGE is a precise and highly sensitive somatic SNV, MNV and small INDEL caller
//

include { SAGE_GERMLINE } from '../../../modules/local/sage/germline/main'
include { SAGE_SOMATIC  } from '../../../modules/local/sage/somatic/main'

workflow SAGE_CALLING {
    take:
    // Sample data
    ch_inputs                    // channel: [mandatory] [ meta ]
    ch_tumor_bam                 // channel: [mandatory] [ meta, bam, bai ]
    ch_normal_bam                // channel: [mandatory] [ meta, bam, bai ]
    ch_donor_bam                 // channel: [mandatory] [ meta, bam, bai ]
    ch_tumor_tsv                 // channel: [mandatory] [ meta, redux_tsv, ... ]
    ch_normal_tsv                // channel: [mandatory] [ meta, redux_tsv, ... ]
    ch_donor_tsv                 // channel: [mandatory] [ meta, redux_tsv, ... ]

    // Reference data
    genome_fasta                 // channel: [mandatory] /path/to/genome_fasta
    genome_version               // channel: [mandatory] genome version
    genome_fai                   // channel: [mandatory] /path/to/genome_fai
    genome_dict                  // channel: [mandatory] /path/to/genome_dict
    sage_pon                     // channel: [mandatory] /path/to/sage_pon
    sage_known_hotspots_somatic  // channel: [mandatory] /path/to/sage_known_hotspots_somatic
    sage_known_hotspots_germline // channel: [optional]  /path/to/sage_known_hotspots_germline
    sage_highconf_regions        // channel: [mandatory] /path/to/sage_highconf_regions
    segment_mappability          // channel: [mandatory] /path/to/segment_mappability
    driver_gene_panel            // channel: [mandatory] /path/to/driver_gene_panel
    ensembl_data_resources       // channel: [mandatory] /path/to/ensembl_data_resources/
    gnomad_resource              // channel: [mandatory] /path/to/gnomad_resource

    // Params
    sequencing_type              // string:  [mandatory] sequencing type
    enable_germline              // boolean: [mandatory] Enable germline mode
    targeted_mode                // boolean: [mandatory] Set targeted mode

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Sort inputs
    // channel: runnable: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, [redux_tsv, ...] ]
    // channel: skip: [ meta ]

    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_tumor_bam, ch_tumor_tsv,
        ch_normal_bam, ch_normal_tsv,
        ch_donor_bam, ch_donor_tsv,
    )
        .map { meta,
            tumor_bam , tumor_bai , tumor_bqr_tsv , tumor_dup_freq_tsv , tumor_jitter_tsv , tumor_ms_tsv,
            normal_bam, normal_bai, normal_bqr_tsv, normal_dup_freq_tsv, normal_jitter_tsv, normal_ms_tsv,
            donor_bam , donor_bai , donor_bqr_tsv , donor_dup_freq_tsv , donor_jitter_tsv , donor_ms_tsv ->

            tumor_bam = Inputs.preferUserProvidedInput(tumor_bam, meta, Constants.INPUT.BAM_REDUX_DNA_TUMOR)
            tumor_bai = Inputs.preferPipelineOutput(tumor_bai, meta, Constants.INPUT.BAI_DNA_TUMOR)

            normal_bam = Inputs.preferUserProvidedInput(normal_bam, meta, Constants.INPUT.BAM_REDUX_DNA_NORMAL)
            normal_bai = Inputs.preferPipelineOutput(normal_bai, meta, Constants.INPUT.BAI_DNA_NORMAL)

            donor_bam = Inputs.preferUserProvidedInput(donor_bam, meta, Constants.INPUT.BAM_REDUX_DNA_DONOR)
            donor_bai = Inputs.preferPipelineOutput(donor_bai, meta, Constants.INPUT.BAI_DNA_DONOR)

            def redux_tsvs = [
                Inputs.preferPipelineOutput(tumor_bqr_tsv, meta, Constants.INPUT.REDUX_BQR_TSV_TUMOR),
                Inputs.preferPipelineOutput(tumor_jitter_tsv, meta, Constants.INPUT.REDUX_JITTER_TSV_TUMOR),
                Inputs.preferPipelineOutput(tumor_ms_tsv, meta, Constants.INPUT.REDUX_MS_TSV_TUMOR),

                Inputs.preferPipelineOutput(normal_bqr_tsv, meta, Constants.INPUT.REDUX_BQR_TSV_NORMAL),
                Inputs.preferPipelineOutput(normal_jitter_tsv, meta, Constants.INPUT.REDUX_JITTER_TSV_NORMAL),
                Inputs.preferPipelineOutput(normal_ms_tsv, meta, Constants.INPUT.REDUX_MS_TSV_NORMAL),

                Inputs.preferPipelineOutput(donor_bqr_tsv, meta, Constants.INPUT.REDUX_BQR_TSV_DONOR),
                Inputs.preferPipelineOutput(donor_jitter_tsv, meta, Constants.INPUT.REDUX_JITTER_TSV_DONOR),
                Inputs.preferPipelineOutput(donor_ms_tsv, meta, Constants.INPUT.REDUX_MS_TSV_DONOR),
            ]

            redux_tsvs = redux_tsvs.findAll { it != [] }

            return [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, redux_tsvs ]
        }
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, redux_tsvs ->
            runnable: tumor_bam
            skip: true
                return meta
        }

    //
    // MODULE: SAGE germline
    //
    // Select inputs that are eligible to run
    // channel: runnable: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, [redux_tsv, ...] ]
    // channel: skip: [ meta ]
    ch_inputs_germline_sorted = ch_inputs_sorted.runnable
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, redux_tsvs ->
            def has_tumor_normal = tumor_bam && normal_bam
            def has_existing = Inputs.hasExistingInput(meta, Constants.INPUT.SAGE_VCF_NORMAL)

            runnable: has_tumor_normal && !has_existing && enable_germline
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_sage, tumor_bam, normal_bam, tumor_bai, normal_bai, [redux_tsv, ...] ]
    ch_sage_germline_inputs = ch_inputs_germline_sorted.runnable
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, redux_tsvs ->

            def meta_sage = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Inputs.getTumorDnaSampleName(meta),
                normal_id: Inputs.getNormalDnaSampleName(meta),
            ]

            return [meta_sage, tumor_bam, normal_bam, tumor_bai, normal_bai, redux_tsvs]
        }

    // Run process
    SAGE_GERMLINE(
        ch_sage_germline_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        sage_known_hotspots_germline,
        sage_highconf_regions,
        driver_gene_panel,
        ensembl_data_resources,
        sequencing_type,
        targeted_mode,
    )

    ch_versions = ch_versions.mix(SAGE_GERMLINE.out.versions)

    //
    // MODULE: SAGE somatic
    //
    // Select inputs that are eligible to run
    // channel: runnable: [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, [redux_tsv, ...] ]
    // channel: skip: [ meta ]
    ch_inputs_somatic_sorted = ch_inputs_sorted.runnable
        .branch { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, redux_tsvs ->
            def has_tumor = tumor_bam
            def has_existing = Inputs.hasExistingInput(meta, Constants.INPUT.SAGE_VCF_TUMOR)

            runnable: has_tumor && !has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: tumor/normal: [ meta_sage, tumor_bam, normal_bam, donor_bam, tumor_bai, normal_bai, donor_bai, [redux_tsv, ...] ]
    // channel: tumor only: [ meta_sage, tumor_bam, [], tumor_bai, [], [redux_tsv, ...] ]
    ch_sage_somatic_inputs = ch_inputs_somatic_sorted.runnable
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, donor_bam, donor_bai, redux_tsvs ->

            def meta_sage = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Inputs.getTumorDnaSampleName(meta),
                donor_id: Inputs.getDonorDnaSampleName(meta),
            ]

            if (normal_bam) {
                meta_sage.normal_id = Inputs.getNormalDnaSampleName(meta)
            }

            if (donor_bam) {
                meta_sage.donor_id = Inputs.getDonorDnaSampleName(meta)
            }

            return [meta_sage, tumor_bam, normal_bam, donor_bam, tumor_bai, normal_bai, donor_bai, redux_tsvs]
        }

    // Run process
    SAGE_SOMATIC(
        ch_sage_somatic_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        sage_pon,
        sage_known_hotspots_somatic,
        sage_highconf_regions,
        driver_gene_panel,
        ensembl_data_resources,
        gnomad_resource,
        sequencing_type,
        targeted_mode,
    )

    ch_versions = ch_versions.mix(SAGE_SOMATIC.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, sage_vcf, sage_tbi ]
    ch_somatic_vcf_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(SAGE_SOMATIC.out.vcf, ch_inputs),
            PlaceholderChannels.vcfTbi(ch_inputs_somatic_sorted.skip),
            PlaceholderChannels.vcfTbi(ch_inputs_sorted.skip),
        )

    // channel: [ meta, sage_vcf, sage_tbi ]
    ch_germline_vcf_out = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(SAGE_GERMLINE.out.vcf, ch_inputs),
            PlaceholderChannels.vcfTbi(ch_inputs_germline_sorted.skip),
            PlaceholderChannels.vcfTbi(ch_inputs_sorted.skip),
        )

    // channel: [ meta, sage_dir ]
    ch_somatic_dir = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(SAGE_SOMATIC.out.sage_dir, ch_inputs),
            PlaceholderChannels.toolDir(ch_inputs_somatic_sorted.skip),
            PlaceholderChannels.toolDir(ch_inputs_sorted.skip),
        )

    // channel: [ meta, sage_dir ]
    ch_germline_dir = Channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(SAGE_GERMLINE.out.sage_dir, ch_inputs),
            PlaceholderChannels.toolDir(ch_inputs_germline_sorted.skip),
            PlaceholderChannels.toolDir(ch_inputs_sorted.skip),
        )

    emit:
    germline_vcf = ch_germline_vcf_out // channel: [ meta, sage_vcf, sage_tbi ]
    somatic_vcf  = ch_somatic_vcf_out  // channel: [ meta, sage_vcf, sage_tbi ]
    germline_dir = ch_germline_dir     // channel: [ meta, sage_dir ]
    somatic_dir  = ch_somatic_dir      // channel: [ meta, sage_dir ]

    versions     = ch_versions         // channel: [ versions.yml ]
}
