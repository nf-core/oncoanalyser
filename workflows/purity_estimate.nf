import Constants
import Processes
import Utils

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Parse input samplesheet
// NOTE(SW): this is done early and outside of gpars so that we can access synchronously and prior to pipeline execution
inputs = Utils.parseInput(params.input, workflow.stubRun, log)

// Get run config
run_config = WorkflowMain.getRunConfig(params, inputs, log)

// Validate inputs
Utils.validateInput(inputs, run_config, params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
]

// TODO(SW): consider whether we should check for null entries here for errors to be more informative
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { AMBER_PROFILING    } from '../subworkflows/local/amber_profiling'
include { COBALT_PROFILING   } from '../subworkflows/local/cobalt_profiling'
include { PREPARE_REFERENCE  } from '../subworkflows/local/prepare_reference'
include { READ_ALIGNMENT_DNA } from '../subworkflows/local/read_alignment_dna'
include { REDUX_PROCESSING   } from '../subworkflows/local/redux_processing'
include { SAGE_APPEND        } from '../subworkflows/local/sage_append'
include { WISP_ANALYSIS      } from '../subworkflows/local/wisp_analysis'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Get absolute file paths
samplesheet = Utils.getFileObject(params.input)

workflow PURITY_ESTIMATE {
    // Create channel for versions
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Create input channel from parsed CSV
    // channel: [ meta ]
    ch_inputs = Channel.fromList(inputs)

    // Set up reference data, assign more human readable variables
    PREPARE_REFERENCE(
        run_config,
    )
    ref_data = PREPARE_REFERENCE.out
    hmf_data = PREPARE_REFERENCE.out.hmf_data

    ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)

    //
    // SUBWORKFLOW: Run read alignment to generate BAMs
    //
    // channel: [ meta, [bam, ...], [bai, ...] ]
    ch_align_dna_tumor_out = Channel.empty()
    ch_align_dna_normal_out = Channel.empty()
    ch_align_dna_donor_out = Channel.empty()
    ch_align_rna_tumor_out = Channel.empty()
    if (run_config.stages.alignment) {

        READ_ALIGNMENT_DNA(
            ch_inputs,
            ref_data.genome_fasta,
            ref_data.genome_bwamem2_index,
            params.max_fastq_records,
            params.fastp_umi,
            params.fastp_umi_location,
            params.fastp_umi_length,
            params.fastp_umi_skip,
        )

        ch_versions = ch_versions.mix(READ_ALIGNMENT_DNA.out.versions)

        ch_align_dna_tumor_out = ch_align_dna_tumor_out.mix(READ_ALIGNMENT_DNA.out.dna_tumor)
        ch_align_dna_normal_out = ch_align_dna_normal_out.mix(READ_ALIGNMENT_DNA.out.dna_normal)
        ch_align_dna_donor_out = ch_align_dna_donor_out.mix(READ_ALIGNMENT_DNA.out.dna_donor)

    } else {

        ch_align_dna_tumor_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_align_dna_normal_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_align_dna_donor_out = ch_inputs.map { meta -> [meta, [], []] }

    }

    //
    // SUBWORKFLOW: Run REDUX for DNA BAMs
    //
    // channel: [ meta, bam, bai ]
    ch_redux_dna_tumor_out = Channel.empty()
    ch_redux_dna_normal_out = Channel.empty()
    ch_redux_dna_donor_out = Channel.empty()

    // channel: [ meta, dup_freq_tsv, jitter_tsv, ms_tsv, repeat_tsv ]
    ch_redux_dna_tumor_tsv_out = Channel.empty()
    ch_redux_dna_normal_tsv_out = Channel.empty()
    ch_redux_dna_donor_tsv_out = Channel.empty()

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
            params.redux_umi,
            params.redux_umi_duplex_delim,
        )

        ch_versions = ch_versions.mix(REDUX_PROCESSING.out.versions)

        ch_redux_dna_tumor_out = ch_redux_dna_tumor_out.mix(REDUX_PROCESSING.out.dna_tumor)
        ch_redux_dna_normal_out = ch_redux_dna_normal_out.mix(REDUX_PROCESSING.out.dna_normal)
        ch_redux_dna_donor_out = ch_redux_dna_donor_out.mix(REDUX_PROCESSING.out.dna_donor)

        ch_redux_dna_tumor_tsv_out = ch_redux_dna_tumor_tsv_out.mix(REDUX_PROCESSING.out.dna_tumor_tsv)
        ch_redux_dna_normal_tsv_out = ch_redux_dna_normal_tsv_out.mix(REDUX_PROCESSING.out.dna_normal_tsv)
        ch_redux_dna_donor_tsv_out = ch_redux_dna_donor_tsv_out.mix(REDUX_PROCESSING.out.dna_donor_tsv)

    } else {

        ch_redux_dna_tumor_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_redux_dna_normal_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_redux_dna_donor_out = ch_inputs.map { meta -> [meta, [], []] }

        ch_redux_dna_tumor_tsv_out = ch_inputs.map { meta -> [meta, [], [], [], []] }
        ch_redux_dna_normal_tsv_out = ch_inputs.map { meta -> [meta, [], [], [], []] }
        ch_redux_dna_donor_tsv_out = ch_inputs.map { meta -> [meta, [], [], [], []] }

    }

    //
    // SUBWORKFLOW: Run AMBER to obtain b-allele frequencies
    //
    // channel: [ meta, amber_dir ]
    ch_amber_out = Channel.empty()
    if (run_config.stages.amber) {

        AMBER_PROFILING(
            ch_inputs,
            ch_redux_dna_tumor_out,
            ch_redux_dna_normal_out,
            ch_redux_dna_donor_out,
            ref_data.genome_version,
            hmf_data.heterozygous_sites,
            [],  // target_region_bed
        )

        ch_versions = ch_versions.mix(AMBER_PROFILING.out.versions)

        ch_amber_out = ch_amber_out.mix(AMBER_PROFILING.out.amber_dir)

    } else {

        ch_amber_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run COBALT to obtain read ratios
    //
    // channel: [ meta, cobalt_dir ]
    ch_cobalt_out = Channel.empty()
    if (run_config.stages.cobalt) {

        COBALT_PROFILING(
            ch_inputs,
            ch_redux_dna_tumor_out,
            ch_redux_dna_normal_out,
            hmf_data.gc_profile,
            hmf_data.diploid_bed,
            [],  // panel_target_region_normalisation
        )

        ch_versions = ch_versions.mix(COBALT_PROFILING.out.versions)

        ch_cobalt_out = ch_cobalt_out.mix(COBALT_PROFILING.out.cobalt_dir)

    } else {

        ch_cobalt_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Append new sample data to primary SAGE WGS VCF
    //
    // channel: [ meta, sage_append_dir ]
    ch_sage_somatic_append_out = Channel.empty()
    if (run_config.stages.orange) {

        SAGE_APPEND(
            ch_inputs,
            ch_inputs.map { meta -> [meta, []] },  // purple_dir
            ch_redux_dna_tumor_out,
            ch_inputs.map { meta -> [meta, [], []] },  // ch_rna_bam
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            false,  // run_germline
        )

        ch_versions = ch_versions.mix(SAGE_APPEND.out.versions)
        ch_sage_somatic_append_out = ch_sage_somatic_append_out.mix(SAGE_APPEND.out.somatic_dir)

    } else {

        ch_sage_somatic_append_out = ch_inputs.map { meta -> [meta, []] }

    }




    // TODO(SW): rework to accept multiple samples per grouping (patient/subject); currently achievable via multiple groups



    //
    // SUBWORKFLOW: Run WISP to estimate tumor purity
    //
    if (run_config.stages.wisp) {

        WISP_ANALYSIS(
            ch_inputs,
            ch_amber_out,
            ch_cobalt_out,
            ch_sage_somatic_append_out,
            ref_data.genome_fasta,
            ref_data.genome_fai,
        )

        ch_versions = ch_versions.mix(WISP_ANALYSIS.out.versions)

    }

    //
    // TASK: Aggregate software versions
    //
    softwareVersionsToYAML(ch_versions)
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
