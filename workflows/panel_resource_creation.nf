/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.types = true

include { AMBER_PROFILING       } from '../subworkflows/local/amber_profiling'
include { COBALT_NORMALISATION  } from '../subworkflows/local/cobalt_normalisation'
include { COBALT_PROFILING      } from '../subworkflows/local/cobalt_profiling'
include { ISOFOX_NORMALISATION  } from '../subworkflows/local/isofox_normalisation'
include { ISOFOX_QUANTIFICATION } from '../subworkflows/local/isofox_quantification'
include { PAVE_PON_CREATION     } from '../subworkflows/local/pave_pon_creation'
include { PREPARE_OUTPUTS       } from '../subworkflows/local/prepare_outputs'
include { PREPARE_REFERENCE     } from '../subworkflows/local/prepare_reference'
include { READ_ALIGNMENT_DNA    } from '../subworkflows/local/read_alignment_dna'
include { READ_ALIGNMENT_RNA    } from '../subworkflows/local/read_alignment_rna'
include { READ_UMI_PROCESSING   } from '../subworkflows/local/read_umi_processing'
include { REDUX_PROCESSING      } from '../subworkflows/local/redux_processing'
include { SAGE_CALLING          } from '../subworkflows/local/sage_calling'

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { getDnaFastqChannel           } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { getRnaFastqChannel           } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { getPrepConfigFromSamplesheet } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/helpers_parameter'
include { Case                         } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/types_records'
include { getSequencingPlatformPon     } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/utils'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PANEL_RESOURCE_CREATION {
    take:
    inputs: List<Case>
    run_config: Map
    params: Map

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
    def target_regions_bed = params.target_regions_bed != null ? file(params.target_regions_bed) : null
    def driver_gene_panel = params.driver_gene_panel != null ? file(params.driver_gene_panel) : null

    def copy_number_percentiles = params.enable_cn_norm_with_wgs_pct ? hmf_data.map { it.copy_number_percentiles } : null

    def isofox_counts = params.isofox_counts != null ? file(params.isofox_counts) : hmf_data.map { it.isofox_counts }
    def isofox_gene_ids = params.isofox_gene_ids != null ? file(params.isofox_gene_ids) : null
    def isofox_gc_ratios = params.isofox_gc_ratios != null ? file(params.isofox_gc_ratios) : hmf_data.map { it.isofox_gc_ratios }
    def isofox_read_length = params.isofox_read_length != null ? params.isofox_read_length : 93  // targeted default read length

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
        ch_inputs.map { meta -> [meta, null, null] },  // ch_dna_donor
        ref_data.genome_fasta,
        ref_data.genome_version,
        ref_data.genome_fai,
        ref_data.genome_dict,
        hmf_data.map { it.unmap_regions },
        hmf_data.map { it.msi_jitter_sites },
        hmf_data.map { it.msi_model_coefficients },
        hmf_data.map { it.msi_model_error_rates },
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
        hmf_data.map { it.ensembl_data_resources },
        driver_gene_panel,
        hmf_data.map { it.known_fusion_data },
        hmf_data.map { it.isofox_excluded_regions },
        hmf_data.map { it.isofox_gene_distribution },
        hmf_data.map { it.isofox_alt_sj_distribution },
        isofox_counts,
        isofox_gc_ratios,
        null,  // isofox_tpm_norm
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
        ch_inputs.map { meta -> [meta, null] },  // ch_redux_dir_donor
        ref_data.genome_fasta,
        ref_data.genome_version,
        ref_data.genome_fai,
        hmf_data.map { it.heterozygous_sites },
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
        hmf_data.map { it.gc_profile },
        hmf_data.map { it.diploid_bed },
        null,  // panel_target_regions_normalisation
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
        ch_inputs.map { meta -> [meta, null] },  // ch_redux_dir_donor
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

    // channel: [ meta, sage_dir ]
    ch_sage_somatic_dir_out = SAGE_CALLING.out.somatic_dir
    ch_sage_germline_dir_out = SAGE_CALLING.out.germline_dir

    //
    // SUBWORKFLOW: Run COBALT normalisation
    //
    COBALT_NORMALISATION(
        ch_amber_out,
        ch_cobalt_out,
        ref_data.genome_version,
        hmf_data.map { it.gc_profile },
        copy_number_percentiles,
        target_regions_bed,
    )
    ch_cobalt_normalisation_tsv = COBALT_NORMALISATION.out.cobalt_normalisation_tsv

    //
    // SUBWORKFLOW: Run PAVE panel of normals creation
    //
    PAVE_PON_CREATION(
        ch_sage_somatic_dir_out,
        ref_data.genome_version,
    )
    ch_pave_pon_panel_creation_artefacts = PAVE_PON_CREATION.out.pave_pon_panel_creation_artefacts

    //
    // SUBWORKFLOW: Run Isofox TPM normalisation
    //
    ch_isofox_normalisation_csv = channel.empty()
    if (run_config.has_rna) {
        ISOFOX_NORMALISATION(
            ch_isofox_out,
            ref_data.genome_version,
            isofox_gene_ids,
            hmf_data.map { it.isofox_gene_distribution },
        )
        ch_isofox_normalisation_csv = ISOFOX_NORMALISATION.out.isofox_normalisation_csv
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
    PREPARE_OUTPUTS(
        ch_amber_out,
        channel.empty(),  // bamtools_tumor
        channel.empty(),  // bamtools_normal
        channel.empty(),  // chord
        channel.empty(),  // cider
        ch_cobalt_out,
        channel.empty(),  // cuppa
        channel.empty(),  // esvee
        ch_align_rna_tumor_out,
        ch_isofox_out,
        channel.empty(),  // lilac
        channel.empty(),  // linx_germline
        channel.empty(),  // linx_somatic
        channel.empty(),  // linx_somatic_visualiser
        channel.empty(),  // linxreport_html
        channel.empty(),  // multiqc
        channel.empty(),  // neo_annotated_fusions
        channel.empty(),  // neo_finder
        channel.empty(),  // neo_scorer_dir
        channel.empty(),  // orange_json
        channel.empty(),  // orange_pdf
        channel.empty(),  // pave_germline
        channel.empty(),  // pave_somatic
        channel.empty(),  // peach
        channel.empty(),  // purple
        channel.empty(),  // qsee
        ch_redux_tumor_out,
        ch_redux_normal_out,
        channel.empty(),  // redux_donor
        channel.empty(),  // sage_append_somatic
        channel.empty(),  // sage_append_germline
        ch_sage_germline_dir_out,
        ch_sage_somatic_dir_out,
        channel.empty(),  // sage_somatic_visualiser
        channel.empty(),  // sigs
        channel.empty(),  // teal_normal_bam
        channel.empty(),  // teal_tumor_bam
        channel.empty(),  // teal_tsvs
        channel.empty(),  // virusbreakend_tsv
        channel.empty(),  // virusbreakend_vcf
        channel.empty(),  // virusinterpreter
        channel.topic('write_reference_data'),
        channel.topic('command_files'),
        channel.empty(),  // wisp
        ch_cobalt_normalisation_tsv,
        ch_isofox_normalisation_csv,
        ch_pave_pon_panel_creation_artefacts,
    )

    emit:
    results = PREPARE_OUTPUTS.out.results
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
