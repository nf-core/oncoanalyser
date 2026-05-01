/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AMBER_PROFILING       } from '../subworkflows/local/amber_profiling'
include { BAMTOOLS_METRICS      } from '../subworkflows/local/bamtools_metrics'
include { CHORD_PREDICTION      } from '../subworkflows/local/chord_prediction'
include { CIDER_CALLING         } from '../subworkflows/local/cider_calling'
include { COBALT_PROFILING      } from '../subworkflows/local/cobalt_profiling'
include { CUPPA_PREDICTION      } from '../subworkflows/local/cuppa_prediction'
include { ESVEE_CALLING         } from '../subworkflows/local/esvee_calling'
include { ISOFOX_QUANTIFICATION } from '../subworkflows/local/isofox_quantification'
include { LILAC_CALLING         } from '../subworkflows/local/lilac_calling'
include { LINX_ANNOTATION       } from '../subworkflows/local/linx_annotation'
include { LINX_PLOTTING         } from '../subworkflows/local/linx_plotting'
include { NEO_PREDICTION        } from '../subworkflows/local/neo_prediction'
include { ORANGE_REPORTING      } from '../subworkflows/local/orange_reporting'
include { PAVE_ANNOTATION       } from '../subworkflows/local/pave_annotation'
include { PEACH_CALLING         } from '../subworkflows/local/peach_calling'
include { PREPARE_REFERENCE     } from '../subworkflows/local/prepare_reference'
include { PURPLE_CALLING        } from '../subworkflows/local/purple_calling'
include { QSEE_METRICS          } from '../subworkflows/local/qsee_metrics'
include { READ_ALIGNMENT_DNA    } from '../subworkflows/local/read_alignment_dna'
include { READ_ALIGNMENT_RNA    } from '../subworkflows/local/read_alignment_rna'
include { REDUX_PROCESSING      } from '../subworkflows/local/redux_processing'
include { SAGE_APPEND           } from '../subworkflows/local/sage_append'
include { SAGE_CALLING          } from '../subworkflows/local/sage_calling'
include { SAGE_PLOTTING         } from '../subworkflows/local/sage_plotting'
include { SIGS_FITTING          } from '../subworkflows/local/sigs_fitting'
include { TEAL_CHARACTERISATION } from '../subworkflows/local/teal_characterisation'
include { VIRUSBREAKEND_CALLING } from '../subworkflows/local/virusbreakend_calling'

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow WGTS {
    take:
    inputs
    stages

    main:

    // Create channel for versions
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Create input channel from parsed CSV
    // channel: [ meta ]
    ch_inputs = Channel.fromList(inputs)

    // Set up reference data, assign more human readable variables
    PREPARE_REFERENCE(
        false, // prepare_reference_only
        inputs,
        stages,
    )
    ref_data = PREPARE_REFERENCE.out
    hmf_data = PREPARE_REFERENCE.out.hmf_data

    ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)

    // Set GRIDSS config
    gridss_config = params.gridss_config !== null ? file(params.gridss_config) : hmf_data.gridss_config

    //
    // SUBWORKFLOW: Run read alignment to generate BAMs
    //
    // channel: [ meta, [bam, ...], [bai, ...] ]
    ch_align_dna_tumor_out = Channel.empty()
    ch_align_dna_normal_out = Channel.empty()
    ch_align_dna_donor_out = Channel.empty()
    ch_align_rna_tumor_out = Channel.empty()
    if (stages.alignment) {

        READ_ALIGNMENT_DNA(
            ch_inputs,
            ref_data.genome_fasta,
            ref_data.genome_bwamem2_index,
            [],  // known_umis
            params.max_fastq_records,
            false,  // fastp_umi_enabled
            '',     // fastp_umi_location
            0,      // fastp_umi_length
            -1,     // fastp_umi_skip
            false,  // fastq_tools_umi_enabled
            '',     // fastq_tools_umi_duplex_delim
        )

        READ_ALIGNMENT_RNA(
            ch_inputs,
            ref_data.genome_star_index,
            [],    // known_umis
            false, // fastq_tools_umi_enabled
            '',    // fastq_tools_umi_duplex_delim
        )

        ch_versions = ch_versions.mix(
            READ_ALIGNMENT_DNA.out.versions,
            READ_ALIGNMENT_RNA.out.versions,
        )

        ch_align_dna_tumor_out = ch_align_dna_tumor_out.mix(READ_ALIGNMENT_DNA.out.dna_tumor)
        ch_align_dna_normal_out = ch_align_dna_normal_out.mix(READ_ALIGNMENT_DNA.out.dna_normal)
        ch_align_dna_donor_out = ch_align_dna_donor_out.mix(READ_ALIGNMENT_DNA.out.dna_donor)
        ch_align_rna_tumor_out = ch_align_rna_tumor_out.mix(READ_ALIGNMENT_RNA.out.rna_tumor)

    } else {

        ch_align_dna_tumor_out = channels.PlaceholderChannels.bamBai(ch_inputs)
        ch_align_dna_normal_out = channels.PlaceholderChannels.bamBai(ch_inputs)
        ch_align_dna_donor_out = channels.PlaceholderChannels.bamBai(ch_inputs)
        ch_align_rna_tumor_out = channels.PlaceholderChannels.bamBai(ch_inputs)

    }

    //
    // SUBWORKFLOW: Run REDUX for DNA BAMs
    //
    // channel: [ meta, bam, bai ]
    ch_redux_dna_tumor_bam_out = Channel.empty()
    ch_redux_dna_normal_bam_out = Channel.empty()
    ch_redux_dna_donor_bam_out = Channel.empty()

    // channel: [ meta, dir ]
    ch_redux_dna_tumor_dir_out = Channel.empty()
    ch_redux_dna_normal_dir_out = Channel.empty()
    ch_redux_dna_donor_dir_out = Channel.empty()

    if (stages.redux) {

        REDUX_PROCESSING(
            ch_inputs,
            ch_align_dna_tumor_out,
            ch_align_dna_normal_out,
            ch_align_dna_donor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            hmf_data.unmap_regions,
            hmf_data.msi_jitter_sites,
            [],  // msi_model_coefficients
            [],  // msi_model_error_rates
            params.sequencing_type,
            false,  // umi_enable
            '',  // umi_duplex_delim
            false,  // targeted_mode
        )

        ch_versions = ch_versions.mix(REDUX_PROCESSING.out.versions)

        ch_redux_dna_tumor_bam_out = ch_redux_dna_tumor_bam_out.mix(REDUX_PROCESSING.out.dna_tumor_bam)
        ch_redux_dna_normal_bam_out = ch_redux_dna_normal_bam_out.mix(REDUX_PROCESSING.out.dna_normal_bam)
        ch_redux_dna_donor_bam_out = ch_redux_dna_donor_bam_out.mix(REDUX_PROCESSING.out.dna_donor_bam)

        ch_redux_dna_tumor_dir_out = ch_redux_dna_tumor_dir_out.mix(REDUX_PROCESSING.out.dna_tumor_dir)
        ch_redux_dna_normal_dir_out = ch_redux_dna_normal_dir_out.mix(REDUX_PROCESSING.out.dna_normal_dir)
        ch_redux_dna_donor_dir_out = ch_redux_dna_donor_dir_out.mix(REDUX_PROCESSING.out.dna_donor_dir)

    } else {

        ch_redux_dna_tumor_bam_out = channels.PlaceholderChannels.bamBai(ch_inputs)
        ch_redux_dna_normal_bam_out = channels.PlaceholderChannels.bamBai(ch_inputs)
        ch_redux_dna_donor_bam_out = channels.PlaceholderChannels.bamBai(ch_inputs)

        ch_redux_dna_tumor_dir_out = channels.PlaceholderChannels.reduxTsvs(ch_inputs)
        ch_redux_dna_normal_dir_out = channels.PlaceholderChannels.reduxTsvs(ch_inputs)
        ch_redux_dna_donor_dir_out = channels.PlaceholderChannels.reduxTsvs(ch_inputs)

    }

    //
    // SUBWORKFLOW: Run Bam Tools to generate stats required for downstream processes
    //
    // channel: [ meta, metrics_dir ]
    ch_bamtools_somatic_out = Channel.empty()
    ch_bamtools_germline_out = Channel.empty()
    if (stages.bamtools) {

        BAMTOOLS_METRICS(
            ch_inputs,
            ch_redux_dna_tumor_bam_out,
            ch_redux_dna_normal_bam_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            hmf_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
            [], // target_region_bed
        )

        ch_versions = ch_versions.mix(BAMTOOLS_METRICS.out.versions)

        ch_bamtools_somatic_out = ch_bamtools_somatic_out.mix(BAMTOOLS_METRICS.out.somatic)
        ch_bamtools_germline_out = ch_bamtools_germline_out.mix(BAMTOOLS_METRICS.out.germline)

    } else {

        ch_bamtools_somatic_out = channels.PlaceholderChannels.toolDir(ch_inputs)
        ch_bamtools_germline_out = channels.PlaceholderChannels.toolDir(ch_inputs)

    }

    //
    // MODULE: Run Isofox to analyse RNA data
    //

    isofox_counts = params.isofox_counts ? file(params.isofox_counts) : hmf_data.isofox_counts
    isofox_gc_ratios = params.isofox_gc_ratios ? file(params.isofox_gc_ratios) : hmf_data.isofox_gc_ratios
    isofox_read_length = params.isofox_read_length !== null ? params.isofox_read_length : pipeline.Constants.DEFAULT_ISOFOX_READ_LENGTH_WTS

    // channel: [ meta, isofox_dir ]
    ch_isofox_out = Channel.empty()
    if (stages.isofox) {

        ISOFOX_QUANTIFICATION(
            ch_inputs,
            ch_align_rna_tumor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            hmf_data.ensembl_data_resources,
            hmf_data.driver_gene_panel,
            hmf_data.known_fusion_data,
            hmf_data.isofox_gene_distribution,
            hmf_data.isofox_alt_sj_distribution,
            isofox_counts,
            isofox_gc_ratios,
            [],  // isofox_tpm_norm
            params.isofox_functions,
            isofox_read_length,
        )

        ch_versions = ch_versions.mix(ISOFOX_QUANTIFICATION.out.versions)

        ch_isofox_out = ch_isofox_out.mix(ISOFOX_QUANTIFICATION.out.isofox_dir)

    } else {

        ch_isofox_out = channels.PlaceholderChannels.toolDir(ch_inputs)

    }

    //
    // SUBWORKFLOW: Run AMBER to obtain b-allele frequencies
    //
    // channel: [ meta, amber_dir ]
    ch_amber_out = Channel.empty()
    if (stages.amber) {

        AMBER_PROFILING(
            ch_inputs,
            ch_redux_dna_tumor_bam_out,
            ch_redux_dna_normal_bam_out,
            ch_redux_dna_donor_bam_out,
            ref_data.genome_version,
            hmf_data.heterozygous_sites,
            [],  // target_region_bed
            [],  // tumor_min_depth
            false,  // purity_estimate_mode
        )

        ch_versions = ch_versions.mix(AMBER_PROFILING.out.versions)

        ch_amber_out = ch_amber_out.mix(AMBER_PROFILING.out.amber_dir)

    } else {

        ch_amber_out = channels.PlaceholderChannels.toolDir(ch_inputs)

    }

    //
    // SUBWORKFLOW: Run COBALT to obtain read ratios
    //
    // channel: [ meta, cobalt_dir ]
    ch_cobalt_out = Channel.empty()
    if (stages.cobalt) {

        COBALT_PROFILING(
            ch_inputs,
            ch_redux_dna_tumor_bam_out,
            ch_redux_dna_normal_bam_out,
            ref_data.genome_version,
            hmf_data.gc_profile,
            hmf_data.diploid_bed,
            [],  // panel_target_region_normalisation
            false,  // targeted_mode
            false,  // purity_estimate_mode
        )

        ch_versions = ch_versions.mix(COBALT_PROFILING.out.versions)

        ch_cobalt_out = ch_cobalt_out.mix(COBALT_PROFILING.out.cobalt_dir)

    } else {

        ch_cobalt_out = channels.PlaceholderChannels.toolDir(ch_inputs)

    }

    //
    // SUBWORKFLOW: Call structural variants with ESVEE
    //
    // channel: [ meta, esvee_dir ]
    ch_esvee_out = Channel.empty()
    if (stages.esvee) {

        ESVEE_CALLING(
            ch_inputs,
            ch_redux_dna_tumor_bam_out,
            ch_redux_dna_normal_bam_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            ref_data.genome_img,
            hmf_data.known_fusions,
            hmf_data.esvee_pon_breakends,
            hmf_data.esvee_pon_breakpoints,
            hmf_data.decoy_sequences_image,
            hmf_data.repeatmasker_annotations,
            hmf_data.unmap_regions,
            [],  // target_region_bed
            params.sequencing_type,
        )

        ch_versions = ch_versions.mix(ESVEE_CALLING.out.versions)

        ch_esvee_out = ch_esvee_out.mix(ESVEE_CALLING.out.esvee_dir)

    } else {

        ch_esvee_out = channels.PlaceholderChannels.toolDir(ch_inputs)

    }

    //
    // SUBWORKFLOW: Call SNV, MNV, and small INDELS with SAGE
    //
    // channel: [ meta, sage_dir ]
    ch_sage_germline_dir_out = Channel.empty()
    ch_sage_somatic_dir_out = Channel.empty()
    if (stages.sage) {

        SAGE_CALLING(
            ch_inputs,
            ch_redux_dna_tumor_bam_out,
            ch_redux_dna_normal_bam_out,
            ch_redux_dna_donor_bam_out,
            ch_redux_dna_tumor_dir_out,
            ch_redux_dna_normal_dir_out,
            ch_redux_dna_donor_dir_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            hmf_data.sage_pon,
            hmf_data.sage_known_hotspots_somatic,
            hmf_data.sage_known_hotspots_germline,
            hmf_data.sage_highconf_regions,
            hmf_data.segment_mappability,
            hmf_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
            hmf_data.gnomad_resource,
            params.sequencing_type,
            true,  // enable_germline
            false, // targeted_mode
        )

        ch_versions = ch_versions.mix(SAGE_CALLING.out.versions)

        ch_sage_germline_dir_out = ch_sage_germline_dir_out.mix(SAGE_CALLING.out.germline_dir)
        ch_sage_somatic_dir_out = ch_sage_somatic_dir_out.mix(SAGE_CALLING.out.somatic_dir)

    } else {

        ch_sage_germline_dir_out = channels.PlaceholderChannels.toolDir(ch_inputs)
        ch_sage_somatic_dir_out = channels.PlaceholderChannels.toolDir(ch_inputs)

    }

    //
    // SUBWORKFLOW: Annotate variants with PAVE
    //
    // channel: [ meta, pave_dir ]
    ch_pave_germline_out = Channel.empty()
    ch_pave_somatic_out = Channel.empty()
    if (stages.pave) {

        PAVE_ANNOTATION(
            ch_inputs,
            ch_sage_germline_dir_out,
            ch_sage_somatic_dir_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            [],  // sage_pon_artefacts
            hmf_data.sage_pon,
            hmf_data.sage_blocklist_regions,
            hmf_data.sage_blocklist_sites,
            hmf_data.clinvar_annotations,
            hmf_data.segment_mappability,
            hmf_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
            hmf_data.gnomad_resource,
            params.sequencing_type,
        )

        ch_versions = ch_versions.mix(PAVE_ANNOTATION.out.versions)

        ch_pave_germline_out = ch_pave_germline_out.mix(PAVE_ANNOTATION.out.germline)
        ch_pave_somatic_out = ch_pave_somatic_out.mix(PAVE_ANNOTATION.out.somatic)

    } else {

        ch_pave_germline_out = channels.PlaceholderChannels.toolDir(ch_inputs)
        ch_pave_somatic_out = channels.PlaceholderChannels.toolDir(ch_inputs)

    }

    //
    // SUBWORKFLOW: Call CNVs, infer purity and ploidy, and recover low quality SVs with PURPLE
    //
    // channel: [ meta, purple_dir ]
    ch_purple_out = Channel.empty()
    if (stages.purple) {

        PURPLE_CALLING(
            ch_inputs,
            ch_amber_out,
            ch_cobalt_out,
            ch_esvee_out,
            ch_pave_somatic_out,
            ch_pave_germline_out,
            channels.PlaceholderChannels.toolDir(ch_inputs), // redux_dir_tumor
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            hmf_data.gc_profile,
            hmf_data.sage_known_hotspots_somatic,
            hmf_data.sage_known_hotspots_germline,
            hmf_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
            hmf_data.germline_amp_del_freq,
            [],  // target_region_bed
        )

        ch_versions = ch_versions.mix(PURPLE_CALLING.out.versions)

        ch_purple_out = ch_purple_out.mix(PURPLE_CALLING.out.purple_dir)

    } else {

        ch_purple_out = channels.PlaceholderChannels.toolDir(ch_inputs)

    }

    //
    // SUBWORKFLOW: QC metrics
    //
    // channel: [ meta, qsee_dir ]
    ch_qsee_out = Channel.empty()
    if (stages.qsee) {

        QSEE_METRICS(
            ch_inputs,
            ch_redux_dna_tumor_dir_out,
            ch_redux_dna_normal_dir_out,
            ch_bamtools_somatic_out,
            ch_bamtools_germline_out,
            ch_cobalt_out,
            ch_esvee_out,
            ch_purple_out,
            hmf_data.driver_gene_panel,
            hmf_data.qsee_cohort_percentiles,
            params.sequencing_type,
            false,  // targeted_mode
        )

        ch_versions = ch_versions.mix(QSEE_METRICS.out.versions)

        ch_qsee_out = ch_qsee_out.mix(QSEE_METRICS.out.qsee_dir)

    } else {

        ch_qsee_out = channels.PlaceholderChannels.toolDir(ch_inputs)

    }

    //
    // SUBWORKFLOW: Append read data to SAGE VCF
    //
    // channel: [ meta, sage_append_dir ]
    ch_sage_somatic_append_out = Channel.empty()
    ch_sage_germline_append_out = Channel.empty()
    if (stages.sage_append) {

        SAGE_APPEND(
            ch_inputs,
            ch_purple_out,
            channels.PlaceholderChannels.bamBai(ch_inputs),  // ch_tumor_redux_bam
            channels.PlaceholderChannels.reduxTsvs(ch_inputs),  // ch_tumor_redux_tsv
            ch_align_rna_tumor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            params.sequencing_type,
            stages.orange,  // enable_germline [run for ORANGE but not Neo]
            false,  // targeted_mode
        )

        ch_versions = ch_versions.mix(SAGE_APPEND.out.versions)

        ch_sage_somatic_append_out = ch_sage_somatic_append_out.mix(SAGE_APPEND.out.somatic_dir)
        ch_sage_germline_append_out = ch_sage_germline_append_out.mix(SAGE_APPEND.out.germline_dir)

    } else {

        ch_sage_somatic_append_out = channels.PlaceholderChannels.toolDir(ch_inputs)
        ch_sage_germline_append_out = channels.PlaceholderChannels.toolDir(ch_inputs)

    }

    //
    // SUBWORKFLOW: Visualise SAGE variants
    //
    if (stages.sage_vis) {

        SAGE_PLOTTING(
            ch_inputs,
            ch_redux_dna_tumor_bam_out,
            ch_redux_dna_normal_bam_out,
            ch_redux_dna_donor_bam_out,
            ch_redux_dna_tumor_dir_out,
            ch_redux_dna_normal_dir_out,
            ch_redux_dna_donor_dir_out,
            ch_purple_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            hmf_data.sage_pon,
            hmf_data.sage_known_hotspots_somatic,
            hmf_data.sage_highconf_regions,
            hmf_data.ensembl_data_resources,
        )

        ch_versions = ch_versions.mix(SAGE_PLOTTING.out.versions)

    }

    //
    // SUBWORKFLOW: Group structural variants into higher order events with LINX
    //
    // channel: [ meta, linx_annotation_dir ]
    ch_linx_somatic_out = Channel.empty()
    ch_linx_germline_out = Channel.empty()
    if (stages.linx) {

        LINX_ANNOTATION(
            ch_inputs,
            ch_purple_out,
            ref_data.genome_version,
            hmf_data.ensembl_data_resources,
            hmf_data.known_fusion_data,
            hmf_data.driver_gene_panel,
        )

        ch_versions = ch_versions.mix(LINX_ANNOTATION.out.versions)

        ch_linx_somatic_out = ch_linx_somatic_out.mix(LINX_ANNOTATION.out.somatic)
        ch_linx_germline_out = ch_linx_germline_out.mix(LINX_ANNOTATION.out.germline)

    } else {

        ch_linx_somatic_out = channels.PlaceholderChannels.toolDir(ch_inputs)
        ch_linx_germline_out = channels.PlaceholderChannels.toolDir(ch_inputs)

    }

    //
    // SUBWORKFLOW: Visualise LINX annotations
    //
    // channel: [ meta, linx_visualiser_dir ]
    ch_linx_somatic_visualiser_dir_out = Channel.empty()
    if (stages.linx) {

        LINX_PLOTTING(
            ch_inputs,
            ch_linx_somatic_out,
            channels.PlaceholderChannels.toolDir(ch_inputs),  // ch_amber
            channels.PlaceholderChannels.toolDir(ch_inputs),  // ch_cobalt
            channels.PlaceholderChannels.toolDir(ch_inputs),  // ch_purple
            ref_data.genome_version,
            hmf_data.ensembl_data_resources,
        )

        ch_versions = ch_versions.mix(LINX_PLOTTING.out.versions)

        ch_linx_somatic_visualiser_dir_out = ch_linx_somatic_visualiser_dir_out.mix(LINX_PLOTTING.out.visualiser_dir)

    } else {

        ch_linx_somatic_visualiser_dir_out = channels.PlaceholderChannels.toolDir(ch_inputs)

    }

    //
    // SUBWORKFLOW: Run CIDER to identify and annotate CDR3 sequences of IG and TCR loci
    //
    if (stages.cider) {

        CIDER_CALLING(
            ch_inputs,
            ch_redux_dna_tumor_bam_out,
            ch_align_rna_tumor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_dict,
            ref_data.genome_img,
        )

        ch_versions = ch_versions.mix(CIDER_CALLING.out.versions)

    }

    //
    // SUBWORKFLOW: Run Sigs to fit somatic smlv to signature definitions
    //
    // channel: [ meta, sigs_dir ]
    ch_sigs_out = Channel.empty()
    if (stages.sigs) {

        SIGS_FITTING(
            ch_inputs,
            ch_purple_out,
            hmf_data.sigs_signatures,
        )

        ch_versions = ch_versions.mix(SIGS_FITTING.out.versions)

        ch_sigs_out = ch_sigs_out.mix(SIGS_FITTING.out.sigs_dir)

    } else {

        ch_sigs_out = channels.PlaceholderChannels.toolDir(ch_inputs)

    }

    //
    // SUBWORKFLOW: Run CHORD to predict HR deficiency status
    //
    // channel: [ meta, chord_dir ]
    ch_chord_out = Channel.empty()
    if (stages.chord) {

        CHORD_PREDICTION(
            ch_inputs,
            ch_purple_out,
            ref_data.genome_fasta,
            ref_data.genome_fai,
            ref_data.genome_dict,
        )

        ch_versions = ch_versions.mix(CHORD_PREDICTION.out.versions)

        ch_chord_out = ch_chord_out.mix(CHORD_PREDICTION.out.chord_dir)

    } else {

        ch_chord_out = channels.PlaceholderChannels.toolDir(ch_inputs)

    }

    //
    // SUBWORKFLOW: Run LILAC for HLA typing and somatic CNV and SNV calling
    //
    // channel: [ meta, lilac_dir ]
    ch_lilac_out = Channel.empty()
    if (stages.lilac) {

        LILAC_CALLING(
            ch_inputs,
            ch_redux_dna_tumor_bam_out,
            ch_redux_dna_normal_bam_out,
            ch_align_rna_tumor_out,
            ch_purple_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            hmf_data.lilac_resources,
            false,  // targeted_mode,
            params.sequencing_type,
        )

        ch_versions = ch_versions.mix(LILAC_CALLING.out.versions)

        ch_lilac_out = ch_lilac_out.mix(LILAC_CALLING.out.lilac_dir)

    } else {

        ch_lilac_out = channels.PlaceholderChannels.toolDir(ch_inputs)

    }

    //
    // SUBWORKFLOW: Run TEAL for characterisation of telometic regions
    //
    if (stages.teal) {

        TEAL_CHARACTERISATION(
            ch_inputs,
            ch_redux_dna_tumor_bam_out,
            ch_redux_dna_normal_bam_out,
            ch_bamtools_somatic_out,
            ch_bamtools_germline_out,
            ch_cobalt_out,
            ch_purple_out,
            ref_data.genome_version,
            params.sequencing_type,
        )

        ch_versions = ch_versions.mix(TEAL_CHARACTERISATION.out.versions)

    }

    //
    // SUBWORKFLOW: Run VIRUSBreakend and Virus Interpreter to quantify viral content
    //
    // channel: [ meta, virusinterpreter_dir ]
    ch_virusinterpreter_out = Channel.empty()

    // NOTE(LN): Virusbreakend currently broken for SBX and Ultima
    def sequencing_type = pipeline.SequencingType.fromString(params.sequencing_type)
    if (stages.virusinterpreter && sequencing_type == pipeline.SequencingType.ILLUMINA) {

        VIRUSBREAKEND_CALLING(
            ch_inputs,
            ch_redux_dna_tumor_bam_out,
            ch_purple_out,
            ch_bamtools_somatic_out,
            ref_data.genome_fasta,
            ref_data.genome_fai,
            ref_data.genome_dict,
            ref_data.genome_gridss_index,
            hmf_data.virusbreakend_db,
            hmf_data.virus_taxonomy_db,
            hmf_data.virus_reporting_db,
            hmf_data.virus_blocklist_db,
            gridss_config,
        )

        ch_versions = ch_versions.mix(VIRUSBREAKEND_CALLING.out.versions)

        ch_virusinterpreter_out = ch_virusinterpreter_out.mix(VIRUSBREAKEND_CALLING.out.virusinterpreter_dir)

    } else {

        ch_virusinterpreter_out = channels.PlaceholderChannels.toolDir(ch_inputs)

    }

    //
    // SUBWORKFLOW: Run PEACH to call germline haplotypes and report pharmacogenomics
    //
    // channel: [ meta, peach_dir ]
    ch_peach_out = Channel.empty()
    if (stages.peach) {

        PEACH_CALLING(
            ch_inputs,
            ch_purple_out,
            hmf_data.peach_haplotypes,
            hmf_data.peach_haplotype_functions,
            hmf_data.peach_drug_info,
        )

        ch_versions = ch_versions.mix(PEACH_CALLING.out.versions)

        ch_peach_out = ch_peach_out.mix(PEACH_CALLING.out.peach_dir)

    } else {

        ch_peach_out = channels.PlaceholderChannels.toolDir(ch_inputs)

    }

    //
    // SUBWORKFLOW: Run Neo to identify and score neoepitopes
    //
    if (stages.neo) {

        NEO_PREDICTION(
            ch_inputs,
            ch_align_rna_tumor_out,
            ch_isofox_out,
            ch_purple_out,
            ch_sage_somatic_append_out,
            ch_lilac_out,
            ch_linx_somatic_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            hmf_data.ensembl_data_resources,
            hmf_data.neo_resources,
            hmf_data.cohort_tpm_medians,
            isofox_read_length,
        )

        ch_versions = ch_versions.mix(NEO_PREDICTION.out.versions)

    }

    //
    // SUBWORKFLOW: Run CUPPA predict tissue of origin
    //
    // channel: [ meta, cuppa_dir ]
    ch_cuppa_out = Channel.empty()
    if (stages.cuppa) {

        CUPPA_PREDICTION(
            ch_inputs,
            ch_isofox_out,
            ch_purple_out,
            ch_linx_somatic_out,
            ch_virusinterpreter_out,
            ref_data.genome_version,
            hmf_data.cuppa_alt_sj,
            hmf_data.cuppa_classifier,
        )

        ch_versions = ch_versions.mix(CUPPA_PREDICTION.out.versions)

        ch_cuppa_out = ch_cuppa_out.mix(CUPPA_PREDICTION.out.cuppa_dir)

    } else {

        ch_cuppa_out = channels.PlaceholderChannels.toolDir(ch_inputs)

    }

    //
    // SUBWORKFLOW: Run ORANGE to generate static PDF report
    //
    if (stages.orange) {

        ORANGE_REPORTING(
            ch_sage_somatic_dir_out,
            ch_sage_germline_dir_out,
            ch_sage_somatic_append_out,
            ch_sage_germline_append_out,
            ch_purple_out,
            ch_qsee_out,
            ch_linx_somatic_out,
            ch_linx_somatic_visualiser_dir_out,
            ch_linx_germline_out,
            ch_virusinterpreter_out,
            ch_chord_out,
            ch_sigs_out,
            ch_lilac_out,
            ch_cuppa_out,
            ch_peach_out,
            ch_isofox_out,
            ref_data.genome_version,
            hmf_data.disease_ontology,
            false,  // targeted_mode
        )

        ch_versions = ch_versions.mix(ORANGE_REPORTING.out.versions)
    }

    //
    // TASK: Aggregate software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'software_versions.yml',
            sort: true,
            newLine: true,
        )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
