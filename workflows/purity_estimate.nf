/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AMBER_PROFILING                 } from '../subworkflows/local/amber_profiling'
include { COBALT_PROFILING                } from '../subworkflows/local/cobalt_profiling'
include { PREPARE_OUTPUTS_PURITY_ESTIMATE } from '../subworkflows/local/prepare_outputs'
include { PREPARE_REFERENCE               } from '../subworkflows/local/prepare_reference'
include { READ_ALIGNMENT_DNA              } from '../subworkflows/local/read_alignment_dna'
include { READ_UMI_PROCESSING             } from '../subworkflows/local/read_umi_processing'
include { REDUX_PROCESSING                } from '../subworkflows/local/redux_processing'
include { SAGE_APPEND                     } from '../subworkflows/local/sage_append'
include { WISP_ANALYSIS                   } from '../subworkflows/local/wisp_analysis'

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { getDnaFastqChannel } from '../subworkflows/local/utils_nfcore_oncoanalyser_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PURITY_ESTIMATE {
    take:
    inputs
    run_config
    params

    main:
    // Create input channel from parsed CSV
    // channel: [ meta ]
    ch_inputs = channel.fromList(inputs)

    // Get run mode of purity estimate mode
    def purity_estimate_run_mode = Utils.getEnumFromString(params.purity_estimate_mode, Constants.RunMode)
    def targeted_mode = purity_estimate_run_mode == Constants.RunMode.TARGETED
    def wgts_mode = purity_estimate_run_mode == Constants.RunMode.WGTS

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
                ch_inputs.map { meta -> [meta, [:], [], []] },  // ch_rna_fastq
                panel_data.known_umis,
                params.fastp_umi_enabled,
                params.fastp_umi_location,
                params.fastp_umi_length,
                params.fastp_umi_skip,
                false,  // fastq_tools_umi_enabled
                '',  // fastq_tools_umi_delim
            )

            ch_align_dna_input = ch_align_dna_input.mix(READ_UMI_PROCESSING.out.fastq_dna)

        } else {

            ch_align_dna_input = ch_fastq_dna

        }

        READ_ALIGNMENT_DNA(
            ch_inputs,
            ch_align_dna_input,
            ref_data.genome_fasta,
            ref_data.genome_bwamem2_index,
            params.max_fastq_records,
        )

        ch_align_dna_tumor_out = ch_align_dna_tumor_out.mix(READ_ALIGNMENT_DNA.out.tumor)
        ch_align_dna_normal_out = ch_align_dna_normal_out.mix(READ_ALIGNMENT_DNA.out.normal)
        ch_align_dna_donor_out = ch_align_dna_donor_out.mix(READ_ALIGNMENT_DNA.out.donor)

    } else {

        ch_align_dna_tumor_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_align_dna_normal_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_align_dna_donor_out = ch_inputs.map { meta -> [meta, [], []] }

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
            // NOTE(LN): panel specific MSI predictions not used as indels are unimportant for WISP
            [],  // msi_model_coefficients
            [],  // msi_model_error_rates
            params.sequencing_platform,
            targeted_mode,
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
            hmf_data.heterozygous_sites,
            [],  // target_regions_bed
            1,  // tumor_min_depth
            params.sequencing_platform,
            true,  // purity_estimate_mode
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
    if (run_config.stages.cobalt && wgts_mode) {

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
            targeted_mode,
            true,  // purity_estimate_mode
        )

        ch_cobalt_out = ch_cobalt_out.mix(COBALT_PROFILING.out.cobalt_dir)

    } else {

        ch_cobalt_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Append new sample data to primary SAGE WGS VCF
    //
    // channel: [ meta, sage_append_dir ]
    ch_sage_somatic_append_out = channel.empty()
    if (run_config.stages.sage_append) {

        SAGE_APPEND(
            ch_inputs,
            ch_inputs.map { meta -> [meta, []] },  // ch_purple_dir
            ch_redux_tumor_out,
            ch_inputs.map { meta -> [meta, [], []] },  // ch_tumor_rna_aln
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            params.sequencing_platform,
            false,  // enable_germline
            targeted_mode,
        )

        ch_sage_somatic_append_out = ch_sage_somatic_append_out.mix(SAGE_APPEND.out.somatic_dir)

    } else {

        ch_sage_somatic_append_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run WISP to estimate tumor purity
    //
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
    PREPARE_OUTPUTS_PURITY_ESTIMATE()

    emit:
    results = PREPARE_OUTPUTS_PURITY_ESTIMATE.out.results
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
