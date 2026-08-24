/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.types = true

include { AMBER_PROFILING                  } from '../subworkflows/local/amber_profiling'
include { COBALT_PROFILING                 } from '../subworkflows/local/cobalt_profiling'
include { PREPARE_OUTPUTS                } from '../subworkflows/local/prepare_outputs'
include { PREPARE_REFERENCE                } from '../subworkflows/local/prepare_reference'
include { READ_ALIGNMENT_DNA               } from '../subworkflows/local/read_alignment_dna'
include { READ_UMI_PROCESSING              } from '../subworkflows/local/read_umi_processing'
include { REDUX_PROCESSING                 } from '../subworkflows/local/redux_processing'
include { SAGE_APPEND                      } from '../subworkflows/local/sage_append'
include { WISP_ANALYSIS                    } from '../subworkflows/local/wisp_analysis'

include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { RunMode                       } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/types_enums'
include { Case                          } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/types_records'

include { getDnaFastqChannel            } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/helpers_channel'
include { getEnumFromString             } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/utils'
include { getPrepConfigFromSamplesheet  } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline/helpers_parameter'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PURITY_ESTIMATE {
    take:
    inputs: List<Case>
    run_config: Map
    params: Map

    main:
    // Create input channel from parsed CSV
    // channel: [ meta ]
    ch_inputs = channel.fromList(inputs)

    // Get run mode of purity estimate mode
    def purity_estimate_run_mode = getEnumFromString(params.purity_estimate_mode, RunMode)
    def targeted_mode = purity_estimate_run_mode == RunMode.TARGETED
    def wgts_mode = purity_estimate_run_mode == RunMode.WGTS

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

    //
    // SUBWORKFLOW: Run read alignment to generate BAMs
    //
    // channel: [ meta, [aln, ...], [idx, ...] ]
    ch_align_dna_tumor_out = channel.empty()
    ch_align_dna_normal_out = channel.empty()
    ch_align_dna_donor_out = channel.empty()
    if (run_config.stages.alignment) {


        // NOTE(LN): For now we won't support purity estimate mode for panel MSK (i.e. UMI processing with fastq-tools)


        // channel: [ meta, fastq_info, fastq_fwd, fastq_rev ]
        ch_fastq_dna = getDnaFastqChannel(ch_inputs)

        // channel: [ meta, fastq_info, fastq_fwd, fastq_rev ]
        ch_align_dna_input = channel.empty()
        if (params.fastp_umi_enabled || params.fastq_tools_umi_enabled) {

            READ_UMI_PROCESSING(
                ch_inputs,
                ch_fastq_dna,
                ch_inputs.map { meta -> [meta, [:], null, null] },  // ch_rna_fastq
                panel_data.map { it.known_umis },
                params.fastp_umi_enabled,
                params.fastp_umi_location,
                params.fastp_umi_length,
                params.fastp_umi_skip,
                false,  // fastq_tools_umi_enabled
                '',  // fastq_tools_umi_delim
            )

            ch_align_dna_input = READ_UMI_PROCESSING.out.fastq_dna
        } else {

            ch_align_dna_input = ch_fastq_dna

        }

        READ_ALIGNMENT_DNA(
            ch_inputs,
            ch_align_dna_input,
            ref_data.genome_fasta,
            ref_data.genome_bwamem2_index,
            params.max_fastq_records,
            params.sequencing_platform,
        )

        ch_align_dna_tumor_out = READ_ALIGNMENT_DNA.out.tumor
        ch_align_dna_normal_out = READ_ALIGNMENT_DNA.out.normal
        ch_align_dna_donor_out = READ_ALIGNMENT_DNA.out.donor
    } else {

        ch_align_dna_tumor_out = ch_inputs.map { meta -> [meta, null, null] }
        ch_align_dna_normal_out = ch_inputs.map { meta -> [meta, null, null] }
        ch_align_dna_donor_out = ch_inputs.map { meta -> [meta, null, null] }

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
            // NOTE(LN): panel specific MSI predictions not used as indels are unimportant for WISP
            null,  // msi_model_coefficients
            null,  // msi_model_error_rates
            params.sequencing_platform,
            targeted_mode,
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
    // SUBWORKFLOW: Run AMBER to obtain b-allele frequencies
    //
    // channel: [ meta, amber_dir ]
    ch_amber_out = channel.empty()
    if (run_config.stages.amber && wgts_mode) {

        AMBER_PROFILING(
            ch_inputs,
            ch_redux_tumor_out,
            ch_redux_normal_out,
            ch_redux_donor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            hmf_data.map { it.heterozygous_sites },
            null,  // target_regions_bed
            1,  // tumor_min_depth
            params.sequencing_platform,
            true,  // purity_estimate_mode
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
    if (run_config.stages.cobalt && wgts_mode) {

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
            targeted_mode,
            true,  // purity_estimate_mode
        )

        ch_cobalt_out = COBALT_PROFILING.out.cobalt_dir
    } else {

        ch_cobalt_out = ch_inputs.map { meta -> [meta, null] }

    }

    //
    // SUBWORKFLOW: Append new sample data to primary SAGE WGS VCF
    //
    // channel: [ meta, sage_append_dir ]
    ch_sage_somatic_append_out = channel.empty()
    if (run_config.stages.sage_append) {

        SAGE_APPEND(
            ch_inputs,
            ch_inputs.map { meta -> [meta, null] },  // ch_purple_dir
            ch_redux_tumor_out,
            ch_inputs.map { meta -> [meta, null, null] },  // ch_tumor_rna_aln
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            params.sequencing_platform,
            false,  // enable_germline
            targeted_mode,
            true,  // purity_estimate_mode
        )

        ch_sage_somatic_append_out = SAGE_APPEND.out.somatic_dir
    } else {

        ch_sage_somatic_append_out = ch_inputs.map { meta -> [meta, null] }

    }

    //
    // SUBWORKFLOW: Run WISP to estimate tumor purity
    //
    // channel: [ meta, wisp_dir ]
    ch_wisp_out = channel.empty()
    if (run_config.stages.wisp) {

        WISP_ANALYSIS(
            ch_inputs,
            ch_redux_tumor_out,
            ch_amber_out,
            ch_cobalt_out,
            ch_sage_somatic_append_out,
            ref_data.genome_fasta,
            ref_data.genome_fai,
            targeted_mode,
        )

        ch_wisp_out = WISP_ANALYSIS.out.wisp_dir
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
        channel.empty(),  // align_rna_tumor
        channel.empty(),  // isofox
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
        ch_redux_donor_out,
        ch_sage_somatic_append_out,
        channel.empty(),  // sage_append_germline
        channel.empty(),  // sage_germline
        channel.empty(),  // sage_somatic
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
        ch_wisp_out,
        channel.empty(),  // cobalt_normalisation_tsv
        channel.empty(),  // isofox_normalisation_csv
        channel.empty(),  // pave_pon_panel_creation_artefacts
    )

    emit:
    results = PREPARE_OUTPUTS.out.results
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
