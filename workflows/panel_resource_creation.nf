/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AMBER_PROFILING                         } from '../subworkflows/local/amber_profiling'
include { COBALT_NORMALISATION                    } from '../subworkflows/local/cobalt_normalisation'
include { COBALT_PROFILING                        } from '../subworkflows/local/cobalt_profiling'
include { ISOFOX_NORMALISATION                    } from '../subworkflows/local/isofox_normalisation'
include { ISOFOX_QUANTIFICATION                   } from '../subworkflows/local/isofox_quantification'
include { PAVE_PON_CREATION                       } from '../subworkflows/local/pave_pon_creation'
include { PREPARE_OUTPUTS_PANEL_RESOURCE_CREATION } from '../subworkflows/local/prepare_outputs'
include { PREPARE_REFERENCE                       } from '../subworkflows/local/prepare_reference'
include { READ_ALIGNMENT_DNA                      } from '../subworkflows/local/read_alignment_dna'
include { READ_ALIGNMENT_RNA                      } from '../subworkflows/local/read_alignment_rna'
include { READ_UMI_PROCESSING                     } from '../subworkflows/local/read_umi_processing'
include { REDUX_PROCESSING                        } from '../subworkflows/local/redux_processing'
include { SAGE_CALLING                            } from '../subworkflows/local/sage_calling'

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { getDnaFastqChannel } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline'
include { getRnaFastqChannel } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PANEL_RESOURCE_CREATION {
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
        params.isofox_gene_ids,
        params.target_regions_bed,
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
    def target_regions_bed = params.target_regions_bed != null ? file(params.target_regions_bed) : []
    def driver_gene_panel = params.driver_gene_panel != null ? file(params.driver_gene_panel) : []

    def copy_number_percentiles = params.enable_cn_norm_with_wgs_pct ? hmf_data.copy_number_percentiles : []

    def isofox_counts = params.isofox_counts != null ? file(params.isofox_counts) : hmf_data.isofox_counts
    def isofox_gene_ids = params.isofox_gene_ids != null ? file(params.isofox_gene_ids) : []
    def isofox_gc_ratios = params.isofox_gc_ratios != null ? file(params.isofox_gc_ratios) : hmf_data.isofox_gc_ratios
    def isofox_read_length = params.isofox_read_length != null ? params.isofox_read_length : Constants.DEFAULT_ISOFOX_READ_LENGTH_TARGETED

    //
    // SUBWORKFLOW: Run read alignment to generate BAMs
    //

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

    // channel: [ meta, [aln, ...], [idx, ...] ]
    ch_align_dna_tumor_out = READ_ALIGNMENT_DNA.out.tumor
    ch_align_dna_normal_out = READ_ALIGNMENT_DNA.out.normal
    ch_align_rna_tumor_out = READ_ALIGNMENT_RNA.out.tumor

    //
    // SUBWORKFLOW: Run REDUX for DNA alignments
    //
    REDUX_PROCESSING(
        ch_inputs,
        ch_align_dna_tumor_out,
        ch_align_dna_normal_out,
        ch_inputs.map { meta -> [meta, [], []] },  // ch_dna_donor
        ref_data.genome_fasta,
        ref_data.genome_version,
        ref_data.genome_fai,
        ref_data.genome_dict,
        hmf_data.unmap_regions,
        hmf_data.msi_jitter_sites,
        hmf_data.msi_model_coefficients,
        hmf_data.msi_model_error_rates,
        params.sequencing_platform,
        true,  // targeted_mode
        params.redux_umi_enabled,
        params.redux_umi_duplex_delim,
    )

    // channel: [ meta, redux_dir ]
    ch_redux_tumor_out = REDUX_PROCESSING.out.tumor_dir
    ch_redux_normal_out = REDUX_PROCESSING.out.normal_dir

    //
    // MODULE: Run Isofox to analyse RNA data
    //
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
        [],  // isofox_tpm_norm
        'TRANSCRIPT_COUNTS',
        isofox_read_length,
    )

    // channel: [ meta, isofox_dir ]
    ch_isofox_out = ISOFOX_QUANTIFICATION.out.isofox_dir

    //
    // SUBWORKFLOW: Run AMBER to obtain b-allele frequencies
    //
    AMBER_PROFILING(
        ch_inputs,
        ch_redux_tumor_out,
        ch_redux_normal_out,
        ch_inputs.map { meta -> [meta, []] },  // ch_redux_dir_donor
        ref_data.genome_fasta,
        ref_data.genome_version,
        ref_data.genome_fai,
        hmf_data.heterozygous_sites,
        target_regions_bed,
        2,  // tumor_min_depth
        params.sequencing_platform,
        false,  // purity_estimate_mode
    )

    // channel: [ meta, amber_dir ]
    ch_amber_out = AMBER_PROFILING.out.amber_dir

    //
    // SUBWORKFLOW: Run COBALT to obtain read ratios
    //
    COBALT_PROFILING(
        ch_inputs,
        ch_redux_tumor_out,
        ch_redux_normal_out,
        ref_data.genome_fasta,
        ref_data.genome_version,
        ref_data.genome_fai,
        hmf_data.gc_profile,
        hmf_data.diploid_bed,
        [],  // panel_target_regions_normalisation
        true,  // targeted_mode
        false,  // purity_estimate_mode
    )

    // channel: [ meta, cobalt_dir ]
    ch_cobalt_out = COBALT_PROFILING.out.cobalt_dir

    //
    // SUBWORKFLOW: call SNV, MNV, and small INDELS with SAGE
    //
    SAGE_CALLING(
        ch_inputs,
        ch_redux_tumor_out,
        ch_redux_normal_out,
        ch_inputs.map { meta -> [meta, []] },  // ch_redux_dir_donor
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

    // channel: [ meta, sage_dir ]
    ch_sage_somatic_dir_out = SAGE_CALLING.out.somatic_dir

    //
    // SUBWORKFLOW: Run COBALT normalisation
    //
    COBALT_NORMALISATION(
        ch_amber_out,
        ch_cobalt_out,
        ref_data.genome_version,
        hmf_data.gc_profile,
        copy_number_percentiles,
        target_regions_bed,
    )

    //
    // SUBWORKFLOW: Run PAVE panel of normals creation
    //
    PAVE_PON_CREATION(
        ch_sage_somatic_dir_out,
        ref_data.genome_version,
    )

    //
    // SUBWORKFLOW: Run Isofox TPM normalisation
    //
    ISOFOX_NORMALISATION(
        ch_isofox_out,
        ref_data.genome_version,
        isofox_gene_ids,
        hmf_data.isofox_gene_distribution,
    )

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

    softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'oncoanalyser_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true,
        )

    //
    // SUBWORKFLOW: Prepare outputs for publishing
    //
    PREPARE_OUTPUTS_PANEL_RESOURCE_CREATION()

    emit:
    results = PREPARE_OUTPUTS_PANEL_RESOURCE_CREATION.out.results
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
