/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.types = true

include { AMBER_PROFILING           } from '../subworkflows/local/amber_profiling'
include { BAMTOOLS_METRICS          } from '../subworkflows/local/bamtools_metrics'
include { CIDER_CALLING             } from '../subworkflows/local/cider_calling'
include { COBALT_PROFILING          } from '../subworkflows/local/cobalt_profiling'
include { ESVEE_CALLING             } from '../subworkflows/local/esvee_calling'
include { ISOFOX_QUANTIFICATION     } from '../subworkflows/local/isofox_quantification'
include { LILAC_CALLING             } from '../subworkflows/local/lilac_calling'
include { LINX_ANNOTATION           } from '../subworkflows/local/linx_annotation'
include { LINX_PLOTTING             } from '../subworkflows/local/linx_plotting'
include { MULTIQC_REPORTING         } from '../subworkflows/local/multiqc_reporting'
include { ORANGE_REPORTING          } from '../subworkflows/local/orange_reporting'
include { PAVE_ANNOTATION           } from '../subworkflows/local/pave_annotation'
include { PEACH_CALLING             } from '../subworkflows/local/peach_calling'
include { get_dir_filepaths       } from '../subworkflows/local/prepare_outputs'
include { get_command_log_filepath } from '../subworkflows/local/prepare_outputs'
include { PREPARE_REFERENCE         } from '../subworkflows/local/prepare_reference'
include { PURPLE_CALLING            } from '../subworkflows/local/purple_calling'
include { QSEE_METRICS              } from '../subworkflows/local/qsee_metrics'
include { READ_ALIGNMENT_DNA        } from '../subworkflows/local/read_alignment_dna'
include { READ_ALIGNMENT_RNA        } from '../subworkflows/local/read_alignment_rna'
include { READ_UMI_PROCESSING       } from '../subworkflows/local/read_umi_processing'
include { REDUX_PROCESSING          } from '../subworkflows/local/redux_processing'
include { SAGE_APPEND               } from '../subworkflows/local/sage_append'
include { SAGE_CALLING              } from '../subworkflows/local/sage_calling'
include { SAGE_PLOTTING             } from '../subworkflows/local/sage_plotting'

include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { getDnaFastqChannel            } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/channels'
include { getRnaFastqChannel            } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/channels'
include { getSequencingPlatformPon     } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/utils'
include { getPrepConfigFromSamplesheet  } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/validate_params'
include { getTumorDnaSampleName         } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/accessors'
include { getNormalDnaSampleName        } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/accessors'
include { getDonorDnaSampleName         } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorRnaSampleName         } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/accessors'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TARGETED {
    take:
    inputs
    run_config
    params

    main:
    // Check input path parameters to see if they exist
    def checkPathParamList = [
        params.driver_gene_panel,
        params.isofox_counts,
        params.isofox_gc_ratios,
        params.isofox_tpm_norm,
    ]

    checkPathParamList.each { param -> if (param) { file(param, checkIfExists: true) } }

    // Create input channel from parsed CSV
    // channel: [ meta ]
    ch_inputs = channel.fromList(inputs)

    // Set up reference data, assign more human readable variables
    def prep_config = getPrepConfigFromSamplesheet(run_config)
    PREPARE_REFERENCE(
        prep_config,
        run_config,
        params,
    )

    def ref_data = PREPARE_REFERENCE.out
    def hmf_data = PREPARE_REFERENCE.out.hmf_data
    def panel_data = PREPARE_REFERENCE.out.panel_data

    // Configure selectable reference data and inputs
    def driver_gene_panel = params.driver_gene_panel != null ? file(params.driver_gene_panel) : panel_data.map { it.driver_gene_panel }
    def msi_model_error_rates = panel_data.map { it.msi_model_error_rates } != null ? panel_data.map { it.msi_model_error_rates } : hmf_data.map { it.msi_model_error_rates }

    def isofox_counts = params.isofox_counts != null ? file(params.isofox_counts) : panel_data.map { it.isofox_counts }
    def isofox_gc_ratios = params.isofox_gc_ratios != null ? file(params.isofox_gc_ratios) : panel_data.map { it.isofox_gc_ratios }
    def isofox_tpm_norm = params.isofox_tpm_norm != null ? file(params.isofox_tpm_norm) : panel_data.map { it.isofox_tpm_norm }
    def isofox_read_length = params.isofox_read_length != null ? params.isofox_read_length : Constants.DEFAULT_ISOFOX_READ_LENGTH_TARGETED

    //
    // SUBWORKFLOW: Run read alignment to generate alignments
    //
    // channel: [ meta, [aln, ...], [idx, ...] ]
    ch_align_dna_tumor_out = channel.empty()
    ch_align_dna_normal_out = channel.empty()
    ch_align_dna_donor_out = channel.empty()
    ch_align_rna_tumor_out = channel.empty()

    // channel: [ meta, star_log, rna_md_metrics ]
    ch_align_rna_qc_tumor_out = channel.empty()

    if (run_config.stages.alignment) {


        // NOTE(SW): fastp can be run twice, multiple passes of the FASTQ in some scenarios, typically not computationally
        // expensive in such situations, so separation between umi / split processing maintained


        // channel: [ meta, fastq_info, fastq_fwd, fastq_rev ]
        ch_fastq_dna = getDnaFastqChannel(ch_inputs)
        ch_fastq_rna = getRnaFastqChannel(ch_inputs)

        // channel: [ meta, fastq_info, fastq_fwd, fastq_rev ]
        ch_align_dna_input = channel.empty()
        ch_align_rna_input = channel.empty()
        if (params.fastp_umi_enabled || params.fastq_tools_umi_enabled) {

            READ_UMI_PROCESSING(
                ch_inputs,
                ch_fastq_dna,
                ch_fastq_rna,
                panel_data.map { it.known_umis },
                params.fastp_umi_enabled,
                params.fastp_umi_location,
                params.fastp_umi_length,
                params.fastp_umi_skip,
                params.fastq_tools_umi_enabled,
                params.fastq_tools_umi_delim,
            )

            ch_align_dna_input = READ_UMI_PROCESSING.out.fastq_dna
            ch_align_rna_input = READ_UMI_PROCESSING.out.fastq_rna
        } else {

            ch_align_dna_input = ch_fastq_dna
            ch_align_rna_input = ch_fastq_rna

        }

        READ_ALIGNMENT_DNA(
            ch_inputs,
            ch_align_dna_input,
            ref_data.genome_fasta,
            ref_data.genome_bwamem2_index,
            params.max_fastq_records,
            params.sequencing_platform,
        )

        READ_ALIGNMENT_RNA(
            ch_inputs,
            ch_align_rna_input,
            ref_data.genome_star_index,
        )

        ch_align_dna_tumor_out = READ_ALIGNMENT_DNA.out.tumor
        ch_align_dna_normal_out = READ_ALIGNMENT_DNA.out.normal
        ch_align_dna_donor_out = READ_ALIGNMENT_DNA.out.donor
        ch_align_rna_tumor_out = READ_ALIGNMENT_RNA.out.tumor
        ch_align_rna_qc_tumor_out = READ_ALIGNMENT_RNA.out.qc_files
    } else {

        ch_align_dna_tumor_out = ch_inputs.map { meta -> [meta, null, null] }
        ch_align_dna_normal_out = ch_inputs.map { meta -> [meta, null, null] }
        ch_align_dna_donor_out = ch_inputs.map { meta -> [meta, null, null] }

        ch_align_rna_tumor_out = ch_inputs.map { meta -> [meta, null, null] }
        ch_align_rna_qc_tumor_out = ch_inputs.map { meta -> [meta, null, null] }

    }

    //
    // SUBWORKFLOW: Run REDUX for DNA alignments
    //
    // channel: [ meta, redux_dir ]
    ch_redux_tumor_out = channel.empty()
    ch_redux_normal_out = channel.empty()
    ch_redux_donor_out = channel.empty()
    if (run_config.stages.redux) {

        REDUX_PROCESSING(
            ch_inputs,
            ch_align_dna_tumor_out,
            ch_align_dna_normal_out,
            ch_align_dna_donor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            hmf_data.map { it.unmap_regions },
            hmf_data.map { it.msi_jitter_sites },
            hmf_data.map { it.msi_model_coefficients },
            msi_model_error_rates,
            params.sequencing_platform,
            true,  // targeted_mode
            params.redux_umi_enabled,
            params.redux_umi_duplex_delim,
        )

        ch_redux_tumor_out = REDUX_PROCESSING.out.tumor_dir
        ch_redux_normal_out = REDUX_PROCESSING.out.normal_dir
        ch_redux_donor_out = REDUX_PROCESSING.out.donor_dir
    } else {

        ch_redux_tumor_out = ch_inputs.map { meta -> [meta, null] }
        ch_redux_normal_out = ch_inputs.map { meta -> [meta, null] }
        ch_redux_donor_out = ch_inputs.map { meta -> [meta, null] }

    }

    //
    // SUBWORKFLOW: Run Bam Tools to generate stats required for downstream processes
    //
    // channel: [ meta, bamtools_dir ]
    ch_bamtools_tumor_out = channel.empty()
    ch_bamtools_normal_out = channel.empty()
    if (run_config.stages.bamtools) {

        BAMTOOLS_METRICS(
            ch_inputs,
            ch_redux_tumor_out,
            ch_redux_normal_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            driver_gene_panel,
            hmf_data.map { it.ensembl_data_resources },
            panel_data.map { it.target_regions_bed },
        )

        ch_bamtools_tumor_out = BAMTOOLS_METRICS.out.tumor_dir
        ch_bamtools_normal_out = BAMTOOLS_METRICS.out.normal_dir
    } else {

        ch_bamtools_tumor_out = ch_inputs.map { meta -> [meta, null] }
        ch_bamtools_normal_out = ch_inputs.map { meta -> [meta, null] }

    }

    //
    // MODULE: Run Isofox to analyse RNA data
    //
    // channel: [ meta, isofox_dir ]
    ch_isofox_out = channel.empty()
    if (run_config.stages.isofox) {

        ISOFOX_QUANTIFICATION(
            ch_inputs,
            ch_align_rna_tumor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            hmf_data.map { it.ensembl_data_resources },
            driver_gene_panel,
            hmf_data.map { it.known_fusion_data },
            hmf_data.map { it.isofox_excluded_regions },
            hmf_data.map { it.isofox_gene_distribution },
            hmf_data.map { it.isofox_alt_sj_distribution },
            isofox_counts,
            isofox_gc_ratios,
            isofox_tpm_norm,
            params.isofox_functions,
            isofox_read_length,
        )

        ch_isofox_out = ISOFOX_QUANTIFICATION.out.isofox_dir
    } else {

        ch_isofox_out = ch_inputs.map { meta -> [meta, null] }

    }

    //
    // SUBWORKFLOW: Run AMBER to obtain b-allele frequencies
    //
    // channel: [ meta, amber_dir ]
    ch_amber_out = channel.empty()
    if (run_config.stages.amber) {

        AMBER_PROFILING(
            ch_inputs,
            ch_redux_tumor_out,
            ch_redux_normal_out,
            ch_redux_donor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            hmf_data.map { it.heterozygous_sites },
            panel_data.map { it.target_regions_bed },
            null,  // tumor_min_depth
            params.sequencing_platform,
            false,  // purity_estimate_mode
        )

        ch_amber_out = AMBER_PROFILING.out.amber_dir
    } else {

        ch_amber_out = ch_inputs.map { meta -> [meta, null] }

    }

    //
    // SUBWORKFLOW: Run COBALT to obtain read ratios
    //
    // channel: [ meta, cobalt_dir ]
    ch_cobalt_out = channel.empty()
    if (run_config.stages.cobalt) {

        COBALT_PROFILING(
            ch_inputs,
            ch_redux_tumor_out,
            ch_redux_normal_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            hmf_data.map { it.gc_profile },
            hmf_data.map { it.diploid_bed },
            panel_data.map { it.target_regions_normalisation },
            true,  // targeted_mode
            false,  // purity_estimate_mode
        )

        ch_cobalt_out = COBALT_PROFILING.out.cobalt_dir
    } else {

        ch_cobalt_out = ch_inputs.map { meta -> [meta, null] }

    }

    //
    // SUBWORKFLOW: Call structural variants with ESVEE
    //
    // channel: [ meta, esvee_dir ]
    ch_esvee_out = channel.empty()
    if (run_config.stages.esvee) {

        ESVEE_CALLING(
            ch_inputs,
            ch_redux_tumor_out,
            ch_redux_normal_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            ref_data.genome_img,
            hmf_data.map { it.known_fusions },
            getSequencingPlatformPon(hmf_data, params.sequencing_platform, 'esvee_breakends', log),
            getSequencingPlatformPon(hmf_data, params.sequencing_platform, 'esvee_breakpoints', log),
            channel.value(null),  // decoy_sequences_image is [] (null) for GRCh38; hmf_data.map of null hangs typed broadcast
            hmf_data.map { it.repeatmasker_annotations },
            hmf_data.map { it.unmap_regions },
            panel_data.map { it.target_regions_bed },
            params.sequencing_platform,
        )

        ch_esvee_out = ESVEE_CALLING.out.esvee_dir
    } else {

        ch_esvee_out = ch_inputs.map { meta -> [meta, null] }

    }

    //
    // SUBWORKFLOW: call SNV, MNV, and small INDELS with SAGE
    //
    // channel: [ meta, sage_dir ]
    ch_sage_germline_out = channel.empty()
    ch_sage_somatic_out = channel.empty()
    if (run_config.stages.sage) {

        SAGE_CALLING(
            ch_inputs,
            ch_redux_tumor_out,
            ch_redux_normal_out,
            ch_redux_donor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            getSequencingPlatformPon(hmf_data, params.sequencing_platform, 'sage', log),
            hmf_data.map { it.sage_known_hotspots_somatic },
            hmf_data.map { it.sage_known_hotspots_germline },
            hmf_data.map { it.sage_highconf_regions },
            hmf_data.map { it.segment_mappability },
            driver_gene_panel,
            hmf_data.map { it.ensembl_data_resources },
            hmf_data.map { it.gnomad_resource },
            params.sequencing_platform,
            true,  // targeted_mode
            true,  // enable_germline
        )

        ch_sage_germline_out = SAGE_CALLING.out.germline_dir
        ch_sage_somatic_out = SAGE_CALLING.out.somatic_dir
    } else {

        ch_sage_germline_out = ch_inputs.map { meta -> [meta, null] }
        ch_sage_somatic_out = ch_inputs.map { meta -> [meta, null] }

    }

    //
    // SUBWORKFLOW: Annotate variants with PAVE
    //
    // channel: [ meta, pave_dir ]
    ch_pave_germline_out = channel.empty()
    ch_pave_somatic_out = channel.empty()
    if (run_config.stages.pave) {

        PAVE_ANNOTATION(
            ch_inputs,
            ch_sage_somatic_out,
            ch_sage_germline_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            panel_data.map { it.pon_artefacts },
            getSequencingPlatformPon(hmf_data, params.sequencing_platform, 'sage', log),
            hmf_data.map { it.sage_blocklist_regions },
            hmf_data.map { it.sage_blocklist_sites },
            hmf_data.map { it.clinvar_annotations },
            hmf_data.map { it.segment_mappability },
            driver_gene_panel,
            hmf_data.map { it.ensembl_data_resources },
            hmf_data.map { it.gnomad_resource },
            params.sequencing_platform,
        )

        ch_pave_germline_out = PAVE_ANNOTATION.out.germline_dir
        ch_pave_somatic_out = PAVE_ANNOTATION.out.somatic_dir
    } else {

        ch_pave_germline_out = ch_inputs.map { meta -> [meta, null] }
        ch_pave_somatic_out = ch_inputs.map { meta -> [meta, null] }

    }

    //
    // SUBWORKFLOW: Call CNVs, infer purity and ploidy, and recover low quality SVs with PURPLE
    //
    // channel: [ meta, purple_dir ]
    ch_purple_out = channel.empty()
    if (run_config.stages.purple) {

        PURPLE_CALLING(
            ch_inputs,
            ch_amber_out,
            ch_cobalt_out,
            ch_esvee_out,
            ch_pave_somatic_out,
            ch_pave_germline_out,
            ch_redux_tumor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            hmf_data.map { it.gc_profile },
            hmf_data.map { it.sage_known_hotspots_somatic },
            hmf_data.map { it.sage_known_hotspots_germline },
            driver_gene_panel,
            hmf_data.map { it.ensembl_data_resources },
            hmf_data.map { it.germline_amp_del_freq },
            panel_data.map { it.target_regions_bed },
        )

        ch_purple_out = PURPLE_CALLING.out.purple_dir
    } else {

        ch_purple_out = ch_inputs.map { meta -> [meta, null] }

    }

    //
    // SUBWORKFLOW: Collect and sumamrise QC metrics with Qsee
    //
    // channel: [ meta, qsee_dir ]
    ch_qsee_out = channel.empty()
    if (run_config.stages.qsee) {

        QSEE_METRICS(
            ch_inputs,
            ch_redux_tumor_out,
            ch_redux_normal_out,
            ch_bamtools_tumor_out,
            ch_bamtools_normal_out,
            ch_cobalt_out,
            ch_esvee_out,
            ch_purple_out,
            driver_gene_panel,
            hmf_data.map { it.qsee_cohort_percentiles },
            params.sequencing_platform,
            true,  // targeted_mode
        )

        ch_qsee_out = QSEE_METRICS.out.qsee_dir
    } else {

        ch_qsee_out = ch_inputs.map { meta -> [meta, null] }

    }

    //
    // SUBWORKFLOW: Append read data to SAGE VCF
    //
    // channel: [ meta, sage_append_dir ]
    ch_sage_append_somatic_out = channel.empty()
    ch_sage_append_germline_out = channel.empty()
    if (run_config.stages.sage_append) {

        SAGE_APPEND(
            ch_inputs,
            ch_purple_out,
            ch_inputs.map { meta -> [meta, null] },  // ch_redux_dir_tumor
            ch_align_rna_tumor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            params.sequencing_platform,
            true,  // enable_germline
            true,  // targeted_mode
        )

        ch_sage_append_somatic_out = SAGE_APPEND.out.somatic_dir
        ch_sage_append_germline_out = SAGE_APPEND.out.germline_dir
    } else {

        ch_sage_append_somatic_out = ch_inputs.map { meta -> [meta, null] }
        ch_sage_append_germline_out = ch_inputs.map { meta -> [meta, null] }

    }

    //
    // SUBWORKFLOW: Visualise SAGE variants
    //
    // channel: [ meta, sage_visualiser_dir ]
    ch_sage_somatic_visualiser_out = channel.empty()
    if (run_config.stages.sage_visualiser) {

        SAGE_PLOTTING(
            ch_inputs,
            ch_redux_tumor_out,
            ch_redux_normal_out,
            ch_redux_donor_out,
            ch_purple_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            getSequencingPlatformPon(hmf_data, params.sequencing_platform, 'sage', log),
            hmf_data.map { it.sage_known_hotspots_somatic },
            hmf_data.map { it.sage_highconf_regions },
            hmf_data.map { it.ensembl_data_resources },
            true,  // targeted_mode
        )

        ch_sage_somatic_visualiser_out = SAGE_PLOTTING.out.visualiser_dir
    } else {

        ch_sage_somatic_visualiser_out = ch_inputs.map { meta -> [meta, null] }

    }

    //
    // SUBWORKFLOW: Group structural variants into higher order events with LINX
    //
    // channel: [ meta, linx_annotation_dir ]
    ch_linx_somatic_out = channel.empty()
    ch_linx_germline_out = channel.empty()
    if (run_config.stages.linx) {

        LINX_ANNOTATION(
            ch_inputs,
            ch_purple_out,
            ref_data.genome_version,
            hmf_data.map { it.ensembl_data_resources },
            hmf_data.map { it.known_fusion_data },
            driver_gene_panel,
        )

        ch_linx_somatic_out = LINX_ANNOTATION.out.somatic_dir
        ch_linx_germline_out = LINX_ANNOTATION.out.germline_dir
    } else {

        ch_linx_somatic_out = ch_inputs.map { meta -> [meta, null] }
        ch_linx_germline_out = ch_inputs.map { meta -> [meta, null] }

    }

    //
    // SUBWORKFLOW: Visualise LINX annotations
    //
    // channel: [ meta, linx_visualiser_dir ]
    ch_linx_somatic_visualiser_out = channel.empty()
    // channel: [ meta, linxreport_html ]
    ch_linxreport_html_out = channel.empty()
    if (run_config.stages.linx) {

        LINX_PLOTTING(
            ch_inputs,
            ch_linx_somatic_out,
            ch_amber_out,
            ch_cobalt_out,
            ch_purple_out,
            ref_data.genome_version,
            hmf_data.map { it.ensembl_data_resources },
        )

        ch_linx_somatic_visualiser_out = LINX_PLOTTING.out.visualiser_dir
        ch_linxreport_html_out = LINX_PLOTTING.out.linxreport_html
    } else {

        ch_linx_somatic_visualiser_out = ch_inputs.map { meta -> [meta, null] }

    }

    //
    // SUBWORKFLOW: Run CIDER to identify and annotate CDR3 sequences of IG and TCR loci
    //
    // channel: [ meta, cider_results ]
    ch_cider_out = channel.empty()
    if (run_config.stages.cider) {

        CIDER_CALLING(
            ch_inputs,
            ch_redux_tumor_out,
            ch_align_rna_tumor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            ref_data.genome_img,
        )

        ch_cider_out = CIDER_CALLING.out.cider_results
    }

    //
    // SUBWORKFLOW: Run LILAC for HLA typing and somatic CNV and SNV calling
    //
    // channel: [ meta, lilac_dir ]
    ch_lilac_out = channel.empty()
    if (run_config.stages.lilac) {

        LILAC_CALLING(
            ch_inputs,
            ch_redux_tumor_out,
            ch_redux_normal_out,
            ch_align_rna_tumor_out,
            ch_purple_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            hmf_data.map { it.lilac_resources },
            params.sequencing_platform,
            true,  // targeted_mode,
        )

        ch_lilac_out = LILAC_CALLING.out.lilac_dir
    } else {

        ch_lilac_out = ch_inputs.map { meta -> [meta, null] }

    }

    //
    // SUBWORKFLOW: Run PEACH to call germline haplotypes and report pharmacogenomics
    //
    // channel: [ meta, peach_dir ]
    ch_peach_out = channel.empty()
    if (run_config.stages.peach) {

        PEACH_CALLING(
            ch_inputs,
            ch_purple_out,
            hmf_data.map { it.peach_haplotypes },
            hmf_data.map { it.peach_haplotype_functions },
            hmf_data.map { it.peach_drug_info },
        )

        ch_peach_out = PEACH_CALLING.out.peach_dir
    } else {

        ch_peach_out = ch_inputs.map { meta -> [meta, null] }

    }

    //
    // SUBWORKFLOW: Run ORANGE to generate static PDF report
    //
    // channel: [ meta, orange_json ] / [ meta, orange_pdf ]
    ch_orange_json_out = channel.empty()
    ch_orange_pdf_out = channel.empty()
    if (run_config.stages.orange) {

        // Create placeholder channels for empty remaining channels
        ch_chord_out = ch_inputs.map { meta -> [meta, null] }
        ch_cuppa_out = ch_inputs.map { meta -> [meta, null] }
        ch_sigs_out = ch_inputs.map { meta -> [meta, null] }
        ch_virusinterpreter_out = ch_inputs.map { meta -> [meta, null] }

        ORANGE_REPORTING(
            ch_sage_somatic_out,
            ch_sage_germline_out,
            ch_sage_append_somatic_out,
            ch_sage_append_germline_out,
            ch_sage_somatic_visualiser_out,
            ch_purple_out,
            ch_qsee_out,
            ch_linx_somatic_out,
            ch_linx_somatic_visualiser_out,
            ch_linx_germline_out,
            ch_virusinterpreter_out,
            ch_chord_out,
            ch_sigs_out,
            ch_lilac_out,
            ch_cuppa_out,
            ch_peach_out,
            ch_isofox_out,
            ref_data.genome_version,
            hmf_data.map { it.disease_ontology },
            params.sequencing_platform,
            true,  // targeted_mode
            params.panel,
        )

        ch_orange_json_out = ORANGE_REPORTING.out.orange_json
        ch_orange_pdf_out = ORANGE_REPORTING.out.orange_pdf
    }

    //
    // TASK: Aggregate software versions
    //
    def topic_versions = channel.topic("versions")
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

    ch_collated_versions = softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'oncoanalyser_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true,
        )

    //
    // SUBWORKFLOW: Run MultiQC to generate HTML report
    //
    // channel: [ multiqc_report ]
    ch_multiqc_out = channel.empty()
    if (run_config.stages.multiqc) {

        MULTIQC_REPORTING(
            ch_bamtools_tumor_out,
            ch_bamtools_normal_out,
            ch_amber_out,
            ch_purple_out,
            ch_align_rna_qc_tumor_out,
            ch_collated_versions,
            params.multiqc_config,
            params.multiqc_methods_description,
            params.multiqc_logo,
        )

        ch_multiqc_out = MULTIQC_REPORTING.out.report
    }

    //
    // SUBWORKFLOW: Prepare outputs for publishing
    //
    ch_results = channel.empty()
        .mix(
            ch_amber_out.flatMap { meta, d ->                     return get_dir_filepaths(meta, d) },
            ch_bamtools_tumor_out.flatMap { meta, d ->            return get_dir_filepaths(meta, d, "bamtools/${getTumorDnaSampleName(meta, primary: true)}") },
            ch_bamtools_normal_out.flatMap { meta, d ->           return get_dir_filepaths(meta, d, "bamtools/${getNormalDnaSampleName(meta)}") },
            ch_cider_out.flatMap { meta, fps ->                   return fps.collect { d -> ["${meta.case_id}/cider/${d.name}", d] } },
            ch_cobalt_out.flatMap { meta, d ->                    return get_dir_filepaths(meta, d) },
            ch_esvee_out.flatMap { meta, d ->                     return get_dir_filepaths(meta, d) },
            ch_align_dna_tumor_out.flatMap { meta, bam, bai ->    return [bam, bai].findAll().collect { d -> ["${meta.case_id}/alignments/${getTumorDnaSampleName(meta, primary: true)}/${d.name}", d] } },
            ch_align_dna_normal_out.flatMap { meta, bam, bai ->   return [bam, bai].findAll().collect { d -> ["${meta.case_id}/alignments/${getNormalDnaSampleName(meta)}/${d.name}", d] } },
            ch_align_dna_donor_out.flatMap { meta, bam, bai ->    return [bam, bai].findAll().collect { d -> ["${meta.case_id}/alignments/${getDonorDnaSampleName(meta)}/${d.name}", d] } },
            ch_align_rna_tumor_out.flatMap { meta, bam, bai ->    return [bam, bai].findAll().collect { d -> ["${meta.case_id}/alignments/${getTumorRnaSampleName(meta)}/${d.name}", d] } },
            ch_isofox_out.flatMap { meta, d ->                    return get_dir_filepaths(meta, d) },
            ch_lilac_out.flatMap { meta, d ->                     return get_dir_filepaths(meta, d) },
            ch_linx_germline_out.flatMap { meta, d ->             return get_dir_filepaths(meta, d, 'linx/germline_annotations') },
            ch_linx_somatic_out.flatMap { meta, d ->              return get_dir_filepaths(meta, d, 'linx/somatic_annotations') },
            ch_linx_somatic_visualiser_out.flatMap { meta, d ->   return get_dir_filepaths(meta, d, 'linx/somatic_plots') },
            ch_linxreport_html_out.map { meta, d ->               return ["${meta.case_id}/linx/${d.name}", d] },
            ch_multiqc_out.flatMap { fps ->                      return fps.collect { d -> [d.name, d] } },
            ch_orange_json_out.filter { meta, d -> d != null }.map { meta, d -> return ["${meta.case_id}/orange/${d.name}", d] },
            ch_orange_pdf_out.filter { meta, d -> d != null }.map { meta, d ->  return ["${meta.case_id}/orange/${d.name}", d] },
            ch_pave_germline_out.flatMap { meta, d ->             return get_dir_filepaths(meta, d, 'pave/germline') },
            ch_pave_somatic_out.flatMap { meta, d ->              return get_dir_filepaths(meta, d, 'pave/somatic') },
            ch_peach_out.flatMap { meta, d ->                     return get_dir_filepaths(meta, d) },
            ch_purple_out.flatMap { meta, d ->                    return get_dir_filepaths(meta, d) },
            ch_qsee_out.flatMap { meta, d ->                      return get_dir_filepaths(meta, d) },
            ch_redux_tumor_out.flatMap { meta, d ->               return get_dir_filepaths(meta, d, "alignments/${getTumorDnaSampleName(meta, primary: true)}") },
            ch_redux_normal_out.flatMap { meta, d ->              return get_dir_filepaths(meta, d, "alignments/${getNormalDnaSampleName(meta)}") },
            ch_redux_donor_out.flatMap { meta, d ->               return get_dir_filepaths(meta, d, "alignments/${getDonorDnaSampleName(meta)}") },
            ch_sage_append_somatic_out.flatMap { meta, d ->       return get_dir_filepaths(meta, d, "sage_append/${getTumorDnaSampleName(meta, primary: true)}") },
            ch_sage_append_germline_out.flatMap { meta, d ->      return get_dir_filepaths(meta, d, "sage_append/${getNormalDnaSampleName(meta)}") },
            ch_sage_germline_out.flatMap { meta, d ->             return get_dir_filepaths(meta, d, 'sage/germline') },
            ch_sage_somatic_out.flatMap { meta, d ->              return get_dir_filepaths(meta, d, 'sage/somatic') },
            ch_sage_somatic_visualiser_out.flatMap { meta, d ->   return get_dir_filepaths(meta, d, 'sage/visualiser') },

            channel.topic('write_reference_data').map { d -> return ["reference_data/${workflow.manifest.version}/", d] },

            channel.topic('command_files').flatMap { f -> get_command_log_filepath(f) }
        )
        .flatMap { meta, d -> return d instanceof Collection ? d.collect { e -> [meta, e] } : [[meta, d]] }

    emit:
    results = ch_results
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
