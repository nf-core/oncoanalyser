import Constants
import Processes
import Utils


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Set defaults and validate parameters
WorkflowOncoanalyser.setParamsDefaults(params, log)
WorkflowOncoanalyser.validateParams(params, log)

// Get parameter summary and also print to console
def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
WorkflowOncoanalyser.paramsSummaryLog(workflow, params, log)

// Set processes to run
processes = Processes.setProcesses(params.mode, log)
processes_include = Processes.getProcessList(params.processes_include, log)
processes_exclude = Processes.getProcessList(params.processes_exclude, log)
Processes.checkIncludeExcludeList(processes_include, processes_exclude, log)

processes.addAll(processes_include)
processes.removeAll(processes_exclude)

run = Constants.Process
    .values()
    .collectEntries { p -> [p.name().toLowerCase(), p in processes] }

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.gridss_config,
    params.linx_gene_id_file,
]

// Conditional requirements
if (run.virusinterpreter) {
    checkPathParamList.add(params.ref_data_virusbreakenddb_path)
}

if (run.lilac && params.ref_data_genome_version == '38' && params.ref_data_genome_type == 'alt') {
    checkPathParamList.add(params.ref_data_hla_slice_bed)
}

// TODO(SW): consider whether we should check for null entries here for errors to be more informative
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Create Path objects for some input files
linx_gene_id_file = params.linx_gene_id_file ? file(params.linx_gene_id_file) : []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOWS
//
include { AMBER_PROFILING       } from '../subworkflows/local/amber_profiling'
include { BAMTOOLS_METRICS      } from '../subworkflows/local/bamtools_metrics'
include { CHORD_PREDICTION      } from '../subworkflows/local/chord_prediction'
include { COBALT_PROFILING      } from '../subworkflows/local/cobalt_profiling'
include { CUPPA_PREDICTION      } from '../subworkflows/local/cuppa_prediction'
include { GRIDSS_CALLING        } from '../subworkflows/local/gridss_calling'
include { GRIDSS_SVPREP_CALLING } from '../subworkflows/local/gridss_svprep_calling'
include { GRIPSS_FILTERING      } from '../subworkflows/local/gripss_filtering'
include { ISOFOX_QUANTIFICATION } from '../subworkflows/local/isofox_quantification'
include { LILAC_CALLING         } from '../subworkflows/local/lilac_calling'
include { LINX_ANNOTATION       } from '../subworkflows/local/linx_annotation'
include { LINX_PLOTTING         } from '../subworkflows/local/linx_plotting'
include { ORANGE_REPORTING      } from '../subworkflows/local/orange_reporting'
include { PAVE_ANNOTATION       } from '../subworkflows/local/pave_annotation'
include { PREPARE_INPUT         } from '../subworkflows/local/prepare_input'
include { PREPARE_REFERENCE     } from '../subworkflows/local/prepare_reference'
include { PURPLE_CALLING        } from '../subworkflows/local/purple_calling'
include { SAGE_CALLING          } from '../subworkflows/local/sage_calling'
include { SIGS_FITTING          } from '../subworkflows/local/sigs_fitting'
include { VIRUSBREAKEND_CALLING } from '../subworkflows/local/virusbreakend_calling'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Get absolute file paths
samplesheet = Utils.getFileObject(params.input)
gridss_config = Utils.getFileObject(params.gridss_config)

workflow WGTS {
    // Create channel for versions
    // channel: [versions.yml]
    ch_versions = Channel.empty()

    // Get inputs from samplesheet, assign more human readable variable
    // channel: [val(meta)]
    PREPARE_INPUT(
        samplesheet,
    )
    ch_inputs = PREPARE_INPUT.out.data

    // Set up reference data, assign more human readable variables
    PREPARE_REFERENCE(
        run,
    )
    ref_data = PREPARE_REFERENCE.out
    hmf_data = PREPARE_REFERENCE.out.hmf_data

    //
    // MODULE: Run Isofox to analyse WTS data
    //
    // channel: [meta, isofox_dir]
    ch_isofox_out = Channel.empty()
    if (run.isofox) {

        ISOFOX_QUANTIFICATION(
            ch_inputs,
            ref_data.genome_fasta,
            ref_data.genome_fai,
            ref_data.genome_version,
            hmf_data.ensembl_data_resources,
            hmf_data.isofox_counts,
            hmf_data.isofox_gc_ratios,
            params.isofox_functions,
        )

        ch_versions = ch_versions.mix(ISOFOX_QUANTIFICATION.out.versions)
        ch_isofox_out = ch_isofox_out.mix(ISOFOX_QUANTIFICATION.out.isofox_dir)
    }

    //
    // SUBWORKFLOW: Run Bam Tools to generate stats required for downstream processes
    //
    // channel: [val(meta), metrics]
    ch_bamtools_somatic_out = Channel.empty()
    ch_bamtools_germline_out = Channel.empty()
    if (run.bamtools) {

        BAMTOOLS_METRICS(
            ch_inputs,
            ref_data.genome_fasta,
            ref_data.genome_version,
            run,
        )

        ch_versions = ch_versions.mix(BAMTOOLS_METRICS.out.versions)
        ch_bamtools_somatic_out = ch_bamtools_somatic_out.mix(BAMTOOLS_METRICS.out.somatic)
        ch_bamtools_germline_out = ch_bamtools_germline_out.mix(BAMTOOLS_METRICS.out.germline)
    }

    //
    // SUBWORKFLOW: Run AMBER to obtain b-allele frequencies
    //
    // channel: [val(meta), amber_dir]
    ch_amber_out = Channel.empty()
    if (run.amber) {

        AMBER_PROFILING(
            ch_inputs,
            ref_data.genome_version,
            hmf_data.heterozygous_sites,
        )

        ch_versions = ch_versions.mix(AMBER_PROFILING.out.versions)
        ch_amber_out = ch_amber_out.mix(AMBER_PROFILING.out.amber_dir)
    }

    //
    // SUBWORKFLOW: Run COBALT to obtain read ratios
    //
    // channel: [val(meta), cobalt_dir]
    ch_cobalt_out = Channel.empty()
    if (run.cobalt) {

        COBALT_PROFILING(
            ch_inputs,
            hmf_data.gc_profile,
        )

        ch_versions = ch_versions.mix(COBALT_PROFILING.out.versions)
        ch_cobalt_out = ch_cobalt_out.mix(COBALT_PROFILING.out.cobalt_dir)
    }

    //
    // SUBWORKFLOW: Call structural variants with GRIDSS
    //
    // channel: [val(meta), gridss_vcf]
    ch_gridss_out = Channel.empty()
    if (run.gridss) {
        if (run.svprep) {

            GRIDSS_SVPREP_CALLING(
                ch_inputs,
                ref_data.genome_fasta,
                ref_data.genome_version,
                ref_data.genome_fai,
                ref_data.genome_dict,
                ref_data.genome_bwa_index,
                ref_data.genome_bwa_index_image,
                ref_data.genome_gridss_index,
                hmf_data.gridss_region_blocklist,
                hmf_data.sv_prep_blocklist,
                hmf_data.known_fusions,
                gridss_config,
            )

            ch_versions = ch_versions.mix(GRIDSS_SVPREP_CALLING.out.versions)
            ch_gridss_out = ch_gridss_out.mix(GRIDSS_SVPREP_CALLING.out.results)

        } else {

            GRIDSS_CALLING(
                ch_inputs,
                ref_data.genome_fasta,
                ref_data.genome_fai,
                ref_data.genome_dict,
                ref_data.genome_bwa_index,
                ref_data.genome_bwa_index_image,
                ref_data.genome_gridss_index,
                hmf_data.gridss_region_blocklist,
                gridss_config,
            )

            ch_versions = ch_versions.mix(GRIDSS_CALLING.out.versions)
            ch_gridss_out = ch_gridss_out.mix(GRIDSS_CALLING.out.results)
        }
    }

    //
    // SUBWORKFLOW: Run GRIPSS to filter GRIDSS SV calls
    //
    // channel: [val(meta), gripss_vcf, gripss_tbi]
    ch_gripss_somatic_out = Channel.empty()
    ch_gripss_germline_out = Channel.empty()
    ch_gripss_somatic_unfiltered_out = Channel.empty()
    if (run.gripss) {

        GRIPSS_FILTERING(
            ch_inputs,
            ch_gridss_out,
            ref_data.genome_fasta,
            ref_data.genome_fai,
            ref_data.genome_version,
            hmf_data.gridss_pon_breakends,
            hmf_data.gridss_pon_breakpoints,
            hmf_data.known_fusions,
            hmf_data.repeatmasker_annotations,
            run,
        )

        ch_versions = ch_versions.mix(GRIPSS_FILTERING.out.versions)
        ch_gripss_somatic_out = ch_gripss_somatic_out.mix(GRIPSS_FILTERING.out.somatic)
        ch_gripss_germline_out = ch_gripss_germline_out.mix(GRIPSS_FILTERING.out.germline)
        ch_gripss_somatic_unfiltered_out = ch_gripss_somatic_unfiltered_out.mix(GRIPSS_FILTERING.out.somatic_unfiltered)
    }

    //
    // SUBWORKFLOW: call SNV, MNV, and small INDELS with SAGE
    //
    // channel: [val(meta), sage_vcf, sage_tbi]
    ch_sage_germline_vcf_out = Channel.empty()
    ch_sage_somatic_vcf_out = Channel.empty()
    // channel: [val(meta), sage_coverage]
    ch_sage_germline_coverage_out = Channel.empty()
    // channel: [val(meta), bqr_plot]
    ch_sage_somatic_tumor_bqr_out = Channel.empty()
    ch_sage_somatic_normal_bqr_out = Channel.empty()
    if (run.sage) {

        SAGE_CALLING(
            ch_inputs,
            ref_data.genome_fasta,
            ref_data.genome_fai,
            ref_data.genome_dict,
            ref_data.genome_version,
            hmf_data.sage_known_hotspots_germline,
            hmf_data.sage_known_hotspots_somatic,
            hmf_data.sage_actionable_panel,
            hmf_data.sage_coverage_panel,
            hmf_data.sage_highconf_regions,
            hmf_data.sage_pon,
            hmf_data.segment_mappability,
            hmf_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(SAGE_CALLING.out.versions)
        ch_sage_germline_vcf_out = ch_sage_germline_vcf_out.mix(SAGE_CALLING.out.germline_vcf)
        ch_sage_germline_coverage_out = ch_sage_germline_coverage_out.mix(SAGE_CALLING.out.germline_coverage)
        ch_sage_somatic_vcf_out = ch_sage_somatic_vcf_out.mix(SAGE_CALLING.out.somatic_vcf)
        ch_sage_somatic_tumor_bqr_out = ch_sage_somatic_tumor_bqr_out.mix(SAGE_CALLING.out.somatic_tumor_bqr)
        ch_sage_somatic_normal_bqr_out = ch_sage_somatic_normal_bqr_out.mix(SAGE_CALLING.out.somatic_normal_bqr)
    }

    //
    // SUBWORKFLOW: Annotate variants with PAVE
    //
    // channel: [val(meta), pave_vcf]
    ch_pave_germline_out = Channel.empty()
    ch_pave_somatic_out = Channel.empty()
    if (run.pave) {

        PAVE_ANNOTATION(
            ch_inputs,
            ch_sage_germline_vcf_out,
            ch_sage_somatic_vcf_out,
            ref_data.genome_fasta,
            ref_data.genome_fai,
            ref_data.genome_version,
            hmf_data.sage_pon,
            hmf_data.sage_blocklist_regions,
            hmf_data.sage_blocklist_sites,
            hmf_data.clinvar_annotations,
            hmf_data.segment_mappability,
            hmf_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
            hmf_data.gnomad_resource,
            run,
        )

        ch_versions = ch_versions.mix(PAVE_ANNOTATION.out.versions)
        ch_pave_germline_out = ch_pave_germline_out.mix(PAVE_ANNOTATION.out.germline)
        ch_pave_somatic_out = ch_pave_somatic_out.mix(PAVE_ANNOTATION.out.somatic)
    }

    //
    // SUBWORKFLOW: Call CNVs, infer purity and ploidy, and recover low quality SVs with PURPLE
    //
    // channel: [val(meta), purple_dir]
    ch_purple_out = Channel.empty()
    if (run.purple) {

        PURPLE_CALLING(
            ch_inputs,
            ch_amber_out,
            ch_cobalt_out,
            ch_pave_somatic_out,
            ch_pave_germline_out,
            ch_gripss_somatic_out,
            ch_gripss_germline_out,
            ch_gripss_somatic_unfiltered_out,
            ref_data.genome_fasta,
            ref_data.genome_fai,
            ref_data.genome_dict,
            ref_data.genome_version,
            hmf_data.gc_profile,
            hmf_data.sage_known_hotspots_somatic,
            hmf_data.sage_known_hotspots_germline,
            hmf_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
            hmf_data.purple_germline_del,
            run,
        )

        ch_versions = ch_versions.mix(PURPLE_CALLING.out.versions)
        ch_purple_out = ch_purple_out.mix(PURPLE_CALLING.out.purple_dir)
    }

    //
    // SUBWORKFLOW: Group structural variants into higher order events with LINX
    //
    // channel: [val(meta), linx_annotation_dir]
    ch_linx_somatic_out = Channel.empty()
    ch_linx_germline_out = Channel.empty()
    // channel: [val(meta), linx_visualiser_dir]
    ch_linx_somatic_plot_out = Channel.empty()
    if (run.linx) {

        LINX_ANNOTATION(
            ch_inputs,
            ch_purple_out,
            ref_data.genome_version,
            hmf_data.ensembl_data_resources,
            hmf_data.known_fusion_data,
            hmf_data.driver_gene_panel,
            linx_gene_id_file,
            run,
        )

        ch_versions = ch_versions.mix(LINX_ANNOTATION.out.versions)
        ch_linx_somatic_out = ch_linx_somatic_out.mix(LINX_ANNOTATION.out.somatic)
        ch_linx_germline_out = ch_linx_germline_out.mix(LINX_ANNOTATION.out.germline)

        LINX_PLOTTING(
            ch_inputs,
            ch_linx_somatic_out,
            ref_data.genome_version,
            hmf_data.ensembl_data_resources,
        )

        ch_versions = ch_versions.mix(LINX_PLOTTING.out.versions)
        ch_linx_somatic_plot_out = ch_linx_somatic_plot_out.mix(LINX_PLOTTING.out.visualiser_dir)
    }

    //
    // SUBWORKFLOW: Run Sigs to fit somatic smlv to signature definitions
    //
    // channel: [val(meta), sigs_dir]
    ch_sigs_out = Channel.empty()
    if (run.sigs) {

        SIGS_FITTING(
            ch_inputs,
            ch_purple_out,
            hmf_data.sigs_signatures,
            run,
        )

        ch_versions = ch_versions.mix(SIGS_FITTING.out.versions)
        ch_sigs_out = ch_sigs_out.mix(SIGS_FITTING.out.sigs_dir)
    }

    //
    // SUBWORKFLOW: Run CHORD to predict HR deficiency status
    //
    // channel: [val(meta), chord_prediction]
    ch_chord_out = Channel.empty()
    if (run.chord) {

        CHORD_PREDICTION(
            ch_inputs,
            ch_purple_out,
            ref_data.genome_version,
            run,
        )

        ch_versions = ch_versions.mix(CHORD_PREDICTION.out.versions)
        ch_chord_out = ch_chord_out.mix(CHORD_PREDICTION.out.prediction)
    }

    //
    // SUBWORKFLOW: Run LILAC for HLA typing and somatic CNV and SNV calling
    //
    // channel: [val(meta), lilac_dir]
    ch_lilac_out = Channel.empty()
    if (run.lilac) {

        // Select HLA slice BED
        if (params.ref_data_genome_type == 'no_alt') {
            ref_data_hla_slice_bed = hmf_data.hla_slice_bed.collect()
        } else if (params.ref_data_genome_type == 'alt' && params.ref_data_genome_version == '38') {
            ref_data_hla_slice_bed = params.ref_data_hla_slice_bed
        } else {
            assert false
        }

        LILAC_CALLING(
            ch_inputs,
            ch_purple_out,
            ref_data.genome_fasta,
            ref_data.genome_fai,
            hmf_data.lilac_resources,
            ref_data_hla_slice_bed,
            run,
        )

        ch_versions = ch_versions.mix(LILAC_CALLING.out.versions)
        ch_lilac_out = ch_lilac_out.mix(LILAC_CALLING.out.lilac_dir)
    }

    //
    // SUBWORKFLOW: Run VIRUSBreakend and Virus Interpreter to quantify viral content
    //
    // channel: [val(meta), virusinterpreter]
    ch_virusinterpreter_out = Channel.empty()
    if (run.virusinterpreter) {

        VIRUSBREAKEND_CALLING(
            ch_inputs,
            ch_purple_out,
            ch_bamtools_somatic_out,
            ref_data.genome_fasta,
            ref_data.genome_fai,
            ref_data.genome_dict,
            ref_data.genome_bwa_index,
            ref_data.genome_bwa_index_image,
            ref_data.genome_gridss_index,
            ref_data.virusbreakenddb,
            hmf_data.virus_taxonomy_db,
            hmf_data.virus_reporting_db,
            run,
            gridss_config,
        )

        ch_versions = ch_versions.mix(VIRUSBREAKEND_CALLING.out.versions)
        ch_virusinterpreter_out = ch_virusinterpreter_out.mix(VIRUSBREAKEND_CALLING.out.virusinterpreter)
    }

    //
    // SUBWORKFLOW: Run CUPPA predict tissue of origin
    //
    // channel: [val(meta), cuppa_dir]
    ch_cuppa_out = Channel.empty()
    if (run.cuppa) {

        CUPPA_PREDICTION(
            ch_inputs,
            ch_isofox_out,
            ch_purple_out,
            ch_linx_somatic_out,
            ch_virusinterpreter_out,
            ref_data.genome_version,
            hmf_data.cuppa_resources,
            run,
        )

        ch_versions = ch_versions.mix(CUPPA_PREDICTION.out.versions)
        ch_cuppa_out = ch_cuppa_out.mix(CUPPA_PREDICTION.out.cuppa_dir)
    }

    //
    // SUBWORKFLOW: Run ORANGE to generate static PDF report
    //
    if (run.orange) {

        ORANGE_REPORTING(
            ch_inputs,
            ch_bamtools_somatic_out,
            ch_bamtools_germline_out,
            ch_sage_somatic_tumor_bqr_out,
            ch_sage_somatic_normal_bqr_out,
            ch_sage_germline_coverage_out,
            ch_purple_out,
            ch_linx_somatic_out,
            ch_linx_somatic_plot_out,
            ch_linx_germline_out,
            ch_virusinterpreter_out,
            ch_chord_out,
            ch_sigs_out,
            ch_lilac_out,
            ch_cuppa_out,
            ch_isofox_out,
            ref_data.genome_version,
            hmf_data.disease_ontology,
            hmf_data.cohort_mapping,
            hmf_data.cohort_percentiles,
            hmf_data.known_fusion_data,
            hmf_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
            hmf_data.alt_sj_distribution,
            hmf_data.gene_exp_distribution,
            run,
        )

        ch_versions = ch_versions.mix(ORANGE_REPORTING.out.versions)
    }

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
