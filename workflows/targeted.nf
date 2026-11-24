/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AMBER_PROFILING          } from '../subworkflows/local/amber_profiling'
include { BAMTOOLS_METRICS         } from '../subworkflows/local/bamtools_metrics'
include { CIDER_CALLING            } from '../subworkflows/local/cider_calling'
include { COBALT_PROFILING         } from '../subworkflows/local/cobalt_profiling'
include { ESVEE_CALLING            } from '../subworkflows/local/esvee_calling'
include { ISOFOX_QUANTIFICATION    } from '../subworkflows/local/isofox_quantification'
include { LILAC_CALLING            } from '../subworkflows/local/lilac_calling'
include { LINX_ANNOTATION          } from '../subworkflows/local/linx_annotation'
include { LINX_PLOTTING            } from '../subworkflows/local/linx_plotting'
include { MULTIQC_REPORTING        } from '../subworkflows/local/multiqc_reporting'
include { ORANGE_REPORTING         } from '../subworkflows/local/orange_reporting'
include { PAVE_ANNOTATION          } from '../subworkflows/local/pave_annotation'
include { PEACH_CALLING            } from '../subworkflows/local/peach_calling'
include { PREPARE_OUTPUTS_TARGETED } from '../subworkflows/local/prepare_outputs'
include { PREPARE_REFERENCE        } from '../subworkflows/local/prepare_reference'
include { PURPLE_CALLING           } from '../subworkflows/local/purple_calling'
include { QSEE_METRICS             } from '../subworkflows/local/qsee_metrics'
include { READ_ALIGNMENT_DNA       } from '../subworkflows/local/read_alignment_dna'
include { READ_ALIGNMENT_RNA       } from '../subworkflows/local/read_alignment_rna'
include { READ_UMI_PROCESSING      } from '../subworkflows/local/read_umi_processing'
include { REDUX_PROCESSING         } from '../subworkflows/local/redux_processing'
include { SAGE_APPEND              } from '../subworkflows/local/sage_append'
include { SAGE_CALLING             } from '../subworkflows/local/sage_calling'
include { SAGE_PLOTTING            } from '../subworkflows/local/sage_plotting'

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { getDnaFastqChannel } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline'
include { getRnaFastqChannel } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline'

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
    def prep_config = WorkflowMain.getPrepConfigFromSamplesheet(run_config)
    PREPARE_REFERENCE(
        prep_config,
        run_config,
        params,
    )

    def ref_data = PREPARE_REFERENCE.out
    def hmf_data = PREPARE_REFERENCE.out.hmf_data
    def panel_data = PREPARE_REFERENCE.out.panel_data

    // Configure selectable reference data and inputs
    def hmf_data_pons = Utils.getSequencingPlatformPons(hmf_data, params.sequencing_platform, log)
    def driver_gene_panel = params.driver_gene_panel != null ? file(params.driver_gene_panel) : panel_data.driver_gene_panel
    def msi_model_error_rates = panel_data.msi_model_error_rates != null ? panel_data.msi_model_error_rates : hmf_data.msi_model_error_rates

    def isofox_counts = params.isofox_counts != null ? file(params.isofox_counts) : panel_data.isofox_counts
    def isofox_gc_ratios = params.isofox_gc_ratios != null ? file(params.isofox_gc_ratios) : panel_data.isofox_gc_ratios
    def isofox_tpm_norm = params.isofox_tpm_norm != null ? file(params.isofox_tpm_norm) : panel_data.isofox_tpm_norm
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
                panel_data.known_umis,
                params.fastp_umi_enabled,
                params.fastp_umi_location,
                params.fastp_umi_length,
                params.fastp_umi_skip,
                params.fastq_tools_umi_enabled,
                params.fastq_tools_umi_delim,
            )

            ch_align_dna_input = ch_align_dna_input.mix(READ_UMI_PROCESSING.out.fastq_dna)
            ch_align_rna_input = ch_align_rna_input.mix(READ_UMI_PROCESSING.out.fastq_rna)

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
        )

        READ_ALIGNMENT_RNA(
            ch_inputs,
            ch_align_rna_input,
            ref_data.genome_star_index,
        )

        ch_align_dna_tumor_out = ch_align_dna_tumor_out.mix(READ_ALIGNMENT_DNA.out.tumor)
        ch_align_dna_normal_out = ch_align_dna_normal_out.mix(READ_ALIGNMENT_DNA.out.normal)
        ch_align_dna_donor_out = ch_align_dna_donor_out.mix(READ_ALIGNMENT_DNA.out.donor)

        ch_align_rna_tumor_out = ch_align_rna_tumor_out.mix(READ_ALIGNMENT_RNA.out.tumor)
        ch_align_rna_qc_tumor_out = ch_align_rna_qc_tumor_out.mix(READ_ALIGNMENT_RNA.out.qc_files)

    } else {

        ch_align_dna_tumor_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_align_dna_normal_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_align_dna_donor_out = ch_inputs.map { meta -> [meta, [], []] }

        ch_align_rna_tumor_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_align_rna_qc_tumor_out = ch_inputs.map { meta -> [meta, [], []] }

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
            hmf_data.unmap_regions,
            hmf_data.msi_jitter_sites,
            hmf_data.msi_model_coefficients,
            msi_model_error_rates,
            params.sequencing_platform,
            true,  // targeted_mode
            params.redux_umi_enabled,
            params.redux_umi_duplex_delim,
        )

        ch_redux_tumor_out = ch_redux_tumor_out.mix(REDUX_PROCESSING.out.tumor_dir)
        ch_redux_normal_out = ch_redux_normal_out.mix(REDUX_PROCESSING.out.normal_dir)
        ch_redux_donor_out = ch_redux_donor_out.mix(REDUX_PROCESSING.out.donor_dir)

    } else {

        ch_redux_tumor_out = ch_inputs.map { meta -> [meta, []] }
        ch_redux_normal_out = ch_inputs.map { meta -> [meta, []] }
        ch_redux_donor_out = ch_inputs.map { meta -> [meta, []] }

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
            hmf_data.ensembl_data_resources,
            panel_data.target_regions_bed,
        )

        ch_bamtools_tumor_out = ch_bamtools_tumor_out.mix(BAMTOOLS_METRICS.out.tumor_dir)
        ch_bamtools_normal_out = ch_bamtools_normal_out.mix(BAMTOOLS_METRICS.out.normal_dir)

    } else {

        ch_bamtools_tumor_out = ch_inputs.map { meta -> [meta, []] }
        ch_bamtools_normal_out = ch_inputs.map { meta -> [meta, []] }

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
            hmf_data.ensembl_data_resources,
            driver_gene_panel,
            hmf_data.known_fusion_data,
            hmf_data.isofox_excluded_regions,
            hmf_data.isofox_gene_distribution,
            hmf_data.isofox_alt_sj_distribution,
            isofox_counts,
            isofox_gc_ratios,
            isofox_tpm_norm,
            params.isofox_functions,
            isofox_read_length,
        )

        ch_isofox_out = ch_isofox_out.mix(ISOFOX_QUANTIFICATION.out.isofox_dir)

    } else {

        ch_isofox_out = ch_inputs.map { meta -> [meta, []] }

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
            hmf_data.heterozygous_sites,
            panel_data.target_regions_bed,
            [],  // tumor_min_depth
            params.sequencing_platform,
            false,  // purity_estimate_mode
        )

        ch_amber_out = ch_amber_out.mix(AMBER_PROFILING.out.amber_dir)

    } else {

        ch_amber_out = ch_inputs.map { meta -> [meta, []] }

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
            hmf_data.gc_profile,
            hmf_data.diploid_bed,
            panel_data.target_regions_normalisation,
            true,  // targeted_mode
            false,  // purity_estimate_mode
        )

        ch_cobalt_out = ch_cobalt_out.mix(COBALT_PROFILING.out.cobalt_dir)

    } else {

        ch_cobalt_out = ch_inputs.map { meta -> [meta, []] }

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
            hmf_data.known_fusions,
            hmf_data_pons.esvee_breakends,
            hmf_data_pons.esvee_breakpoints,
            hmf_data.decoy_sequences_image,
            hmf_data.repeatmasker_annotations,
            hmf_data.unmap_regions,
            panel_data.target_regions_bed,
            params.sequencing_platform,
        )

        ch_esvee_out = ch_esvee_out.mix(ESVEE_CALLING.out.esvee_dir)

    } else {

        ch_esvee_out = ch_inputs.map { meta -> [meta, []] }

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
            hmf_data_pons.sage,
            hmf_data.sage_known_hotspots_somatic,
            hmf_data.sage_known_hotspots_germline,
            hmf_data.sage_highconf_regions,
            hmf_data.segment_mappability,
            driver_gene_panel,
            hmf_data.ensembl_data_resources,
            hmf_data.gnomad_resource,
            params.sequencing_platform,
            true,  // targeted_mode
            true,  // enable_germline
        )

        ch_sage_germline_out = ch_sage_germline_out.mix(SAGE_CALLING.out.germline_dir)
        ch_sage_somatic_out = ch_sage_somatic_out.mix(SAGE_CALLING.out.somatic_dir)

    } else {

        ch_sage_germline_out = ch_inputs.map { meta -> [meta, []] }
        ch_sage_somatic_out = ch_inputs.map { meta -> [meta, []] }

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
            panel_data.pon_artefacts,
            hmf_data_pons.sage,
            hmf_data.sage_blocklist_regions,
            hmf_data.sage_blocklist_sites,
            hmf_data.clinvar_annotations,
            hmf_data.segment_mappability,
            driver_gene_panel,
            hmf_data.ensembl_data_resources,
            hmf_data.gnomad_resource,
            params.sequencing_platform,
        )

        ch_pave_germline_out = ch_pave_germline_out.mix(PAVE_ANNOTATION.out.germline_dir)
        ch_pave_somatic_out = ch_pave_somatic_out.mix(PAVE_ANNOTATION.out.somatic_dir)

    } else {

        ch_pave_germline_out = ch_inputs.map { meta -> [meta, []] }
        ch_pave_somatic_out = ch_inputs.map { meta -> [meta, []] }

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
            hmf_data.gc_profile,
            hmf_data.sage_known_hotspots_somatic,
            hmf_data.sage_known_hotspots_germline,
            driver_gene_panel,
            hmf_data.ensembl_data_resources,
            hmf_data.germline_amp_del_freq,
            panel_data.target_regions_bed,
        )

        ch_purple_out = ch_purple_out.mix(PURPLE_CALLING.out.purple_dir)

    } else {

        ch_purple_out = ch_inputs.map { meta -> [meta, []] }

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
            hmf_data.qsee_cohort_percentiles,
            params.sequencing_platform,
            true,  // targeted_mode
        )

        ch_qsee_out = ch_qsee_out.mix(QSEE_METRICS.out.qsee_dir)

    } else {

        ch_qsee_out = ch_inputs.map { meta -> [meta, []] }

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
            ch_inputs.map { meta -> [meta, []] },  // ch_redux_dir_tumor
            ch_align_rna_tumor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            params.sequencing_platform,
            true,  // enable_germline
            true,  // targeted_mode
        )

        ch_sage_append_somatic_out = ch_sage_append_somatic_out.mix(SAGE_APPEND.out.somatic_dir)
        ch_sage_append_germline_out = ch_sage_append_germline_out.mix(SAGE_APPEND.out.germline_dir)

    } else {

        ch_sage_append_somatic_out = ch_inputs.map { meta -> [meta, []] }
        ch_sage_append_germline_out = ch_inputs.map { meta -> [meta, []] }

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
            hmf_data_pons.sage,
            hmf_data.sage_known_hotspots_somatic,
            hmf_data.sage_highconf_regions,
            hmf_data.ensembl_data_resources,
            true,  // targeted_mode
        )

        ch_sage_somatic_visualiser_out = ch_sage_somatic_visualiser_out.mix(SAGE_PLOTTING.out.visualiser_dir)

    } else {

        ch_sage_somatic_visualiser_out = ch_inputs.map { meta -> [meta, []] }

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
            hmf_data.ensembl_data_resources,
            hmf_data.known_fusion_data,
            driver_gene_panel,
        )

        ch_linx_somatic_out = ch_linx_somatic_out.mix(LINX_ANNOTATION.out.somatic_dir)
        ch_linx_germline_out = ch_linx_germline_out.mix(LINX_ANNOTATION.out.germline_dir)

    } else {

        ch_linx_somatic_out = ch_inputs.map { meta -> [meta, []] }
        ch_linx_germline_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Visualise LINX annotations
    //
    // channel: [ meta, linx_visualiser_dir ]
    ch_linx_somatic_visualiser_out = channel.empty()
    if (run_config.stages.linx) {

        LINX_PLOTTING(
            ch_inputs,
            ch_linx_somatic_out,
            ch_amber_out,
            ch_cobalt_out,
            ch_purple_out,
            ref_data.genome_version,
            hmf_data.ensembl_data_resources,
        )

        ch_linx_somatic_visualiser_out = ch_linx_somatic_visualiser_out.mix(LINX_PLOTTING.out.visualiser_dir)

    } else {

        ch_linx_somatic_visualiser_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run CIDER to identify and annotate CDR3 sequences of IG and TCR loci
    //
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
            hmf_data.lilac_resources,
            params.sequencing_platform,
            true,  // targeted_mode,
        )

        ch_lilac_out = ch_lilac_out.mix(LILAC_CALLING.out.lilac_dir)

    } else {

        ch_lilac_out = ch_inputs.map { meta -> [meta, []] }

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
            hmf_data.peach_haplotypes,
            hmf_data.peach_haplotype_functions,
            hmf_data.peach_drug_info,
        )

        ch_peach_out = ch_peach_out.mix(PEACH_CALLING.out.peach_dir)

    } else {

        ch_peach_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run ORANGE to generate static PDF report
    //
    if (run_config.stages.orange) {

        // Create placeholder channels for empty remaining channels
        ch_chord_out = ch_inputs.map { meta -> [meta, []] }
        ch_cuppa_out = ch_inputs.map { meta -> [meta, []] }
        ch_sigs_out = ch_inputs.map { meta -> [meta, []] }
        ch_virusinterpreter_out = ch_inputs.map { meta -> [meta, []] }

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
            hmf_data.disease_ontology,
            params.sequencing_platform,
            true,  // targeted_mode
            params.panel,
        )

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

        ch_multiqc_out = ch_multiqc_out.mix(MULTIQC_REPORTING.out.report)

    }

    //
    // SUBWORKFLOW: Prepare outputs for publishing
    //
    PREPARE_OUTPUTS_TARGETED()

    emit:
    results = PREPARE_OUTPUTS_TARGETED.out.results
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
