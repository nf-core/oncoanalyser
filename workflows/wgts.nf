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
Utils.validateInput(inputs, run_config, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.isofox_counts,
    params.isofox_gc_ratios,
    params.linx_gene_id_file,
]

// Conditional requirements
if (run_config.stages.gridss) {
    if (params.containsKey('gridss_config')) {
        checkPathParamList.add(params.gridss_config)
    }
}

if (run_config.stages.virusinterpreter) {
    checkPathParamList.add(params.ref_data_virusbreakenddb_path)
}

if (run_config.stages.lilac) {
    if (params.ref_data_genome_version == '38' && params.ref_data_genome_type == 'alt' && params.containsKey('ref_data_hla_slice_bed')) {
        checkPathParamList.add(params.ref_data_hla_slice_bed)
    }
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
include { FLAGSTAT_METRICS      } from '../subworkflows/local/flagstat_metrics'
include { GRIDSS_SVPREP_CALLING } from '../subworkflows/local/gridss_svprep_calling'
include { GRIPSS_FILTERING      } from '../subworkflows/local/gripss_filtering'
include { ISOFOX_QUANTIFICATION } from '../subworkflows/local/isofox_quantification'
include { LILAC_CALLING         } from '../subworkflows/local/lilac_calling'
include { LINX_ANNOTATION       } from '../subworkflows/local/linx_annotation'
include { LINX_PLOTTING         } from '../subworkflows/local/linx_plotting'
include { ORANGE_REPORTING      } from '../subworkflows/local/orange_reporting'
include { PAVE_ANNOTATION       } from '../subworkflows/local/pave_annotation'
include { PREPARE_REFERENCE     } from '../subworkflows/local/prepare_reference'
include { PURPLE_CALLING        } from '../subworkflows/local/purple_calling'
include { SAGE_APPEND           } from '../subworkflows/local/sage_append'
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

workflow WGTS {
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

    // Set GRIDSS config
    gridss_config = params.containsKey('gridss_config') ? file(params.gridss_config) : hmf_data.gridss_config

    //
    // MODULE: Run Isofox to analyse RNA data
    //
    // channel: [ meta, isofox_dir ]
    ch_isofox_out = Channel.empty()
    if (run_config.stages.isofox) {

        isofox_counts = params.isofox_counts ? file(params.isofox_counts) : hmf_data.isofox_counts
        isofox_gc_ratios = params.isofox_gc_ratios ? file(params.isofox_gc_ratios) : hmf_data.isofox_gc_ratios

        ISOFOX_QUANTIFICATION(
            ch_inputs,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            hmf_data.ensembl_data_resources,
            isofox_counts,
            isofox_gc_ratios,
            [],  // isofox_gene_ids
            [],  // isofox_tpm_norm
            params.isofox_functions,
            params.isofox_read_length,
        )

        ch_versions = ch_versions.mix(ISOFOX_QUANTIFICATION.out.versions)

        ch_isofox_out = ch_isofox_out.mix(ISOFOX_QUANTIFICATION.out.isofox_dir)

    } else {

        ch_isofox_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run AMBER to obtain b-allele frequencies
    //
    // channel: [ meta, amber_dir ]
    ch_amber_out = Channel.empty()
    if (run_config.stages.amber) {

        AMBER_PROFILING(
            ch_inputs,
            ref_data.genome_version,
            hmf_data.heterozygous_sites,
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
    // SUBWORKFLOW: Call structural variants with GRIDSS
    //
    // channel: [ meta, gridss_vcf ]
    ch_gridss_out = Channel.empty()
    if (run_config.stages.gridss) {

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

        ch_gridss_out = ch_gridss_out.mix(GRIDSS_SVPREP_CALLING.out.vcf)

    } else {

        ch_gridss_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run GRIPSS to filter GRIDSS SV calls
    //
    // channel: [ meta, gripss_vcf, gripss_tbi ]
    ch_gripss_somatic_out = Channel.empty()
    ch_gripss_germline_out = Channel.empty()
    ch_gripss_somatic_unfiltered_out = Channel.empty()
    if (run_config.stages.gripss) {

        GRIPSS_FILTERING(
            ch_inputs,
            ch_gridss_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            hmf_data.gridss_pon_breakends,
            hmf_data.gridss_pon_breakpoints,
            hmf_data.known_fusions,
            hmf_data.repeatmasker_annotations,
            [],  // target_region_bed
        )

        ch_versions = ch_versions.mix(GRIPSS_FILTERING.out.versions)

        ch_gripss_somatic_out = ch_gripss_somatic_out.mix(GRIPSS_FILTERING.out.somatic)
        ch_gripss_germline_out = ch_gripss_germline_out.mix(GRIPSS_FILTERING.out.germline)
        ch_gripss_somatic_unfiltered_out = ch_gripss_somatic_unfiltered_out.mix(GRIPSS_FILTERING.out.somatic_unfiltered)

    } else {

        ch_gripss_somatic_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_gripss_germline_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_gripss_somatic_unfiltered_out = ch_inputs.map { meta -> [meta, [], []] }

    }

    //
    // SUBWORKFLOW: call SNV, MNV, and small INDELS with SAGE
    //
    // channel: [ meta, sage_vcf, sage_tbi ]
    ch_sage_germline_vcf_out = Channel.empty()
    ch_sage_somatic_vcf_out = Channel.empty()
    // channel: [ meta, sage_dir ]
    ch_sage_germline_dir_out = Channel.empty()
    ch_sage_somatic_dir_out = Channel.empty()
    if (run_config.stages.sage) {

        SAGE_CALLING(
            ch_inputs,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            hmf_data.sage_known_hotspots_germline,
            hmf_data.sage_known_hotspots_somatic,
            hmf_data.sage_actionable_panel,
            hmf_data.sage_coverage_panel,
            hmf_data.sage_highconf_regions,
            hmf_data.segment_mappability,
            hmf_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
        )

        ch_versions = ch_versions.mix(SAGE_CALLING.out.versions)

        ch_sage_germline_vcf_out = ch_sage_germline_vcf_out.mix(SAGE_CALLING.out.germline_vcf)
        ch_sage_somatic_vcf_out = ch_sage_somatic_vcf_out.mix(SAGE_CALLING.out.somatic_vcf)
        ch_sage_germline_dir_out = ch_sage_germline_dir_out.mix(SAGE_CALLING.out.germline_dir)
        ch_sage_somatic_dir_out = ch_sage_somatic_dir_out.mix(SAGE_CALLING.out.somatic_dir)

    } else {

        ch_sage_germline_vcf_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_sage_somatic_vcf_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_sage_germline_dir_out = ch_inputs.map { meta -> [meta, []] }
        ch_sage_somatic_dir_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Annotate variants with PAVE
    //
    // channel: [ meta, pave_vcf ]
    ch_pave_germline_out = Channel.empty()
    ch_pave_somatic_out = Channel.empty()
    if (run_config.stages.pave) {

        PAVE_ANNOTATION(
            ch_inputs,
            ch_sage_germline_vcf_out,
            ch_sage_somatic_vcf_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            hmf_data.sage_pon,
            [],  // sage_pon_artefacts
            hmf_data.sage_blocklist_regions,
            hmf_data.sage_blocklist_sites,
            hmf_data.clinvar_annotations,
            hmf_data.segment_mappability,
            hmf_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
            hmf_data.gnomad_resource,
        )

        ch_versions = ch_versions.mix(PAVE_ANNOTATION.out.versions)

        ch_pave_germline_out = ch_pave_germline_out.mix(PAVE_ANNOTATION.out.germline)
        ch_pave_somatic_out = ch_pave_somatic_out.mix(PAVE_ANNOTATION.out.somatic)

    } else {

        ch_pave_germline_out = ch_inputs.map { meta -> [meta, []] }
        ch_pave_somatic_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Call CNVs, infer purity and ploidy, and recover low quality SVs with PURPLE
    //
    // channel: [ meta, purple_dir ]
    ch_purple_out = Channel.empty()
    if (run_config.stages.purple) {

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
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            hmf_data.gc_profile,
            hmf_data.sage_known_hotspots_somatic,
            hmf_data.sage_known_hotspots_germline,
            hmf_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
            hmf_data.purple_germline_del,
            [],  // target_region_bed
            [],  // target_region_ratios
            [],  // target_region_msi_indels
        )

        ch_versions = ch_versions.mix(PURPLE_CALLING.out.versions)

        ch_purple_out = ch_purple_out.mix(PURPLE_CALLING.out.purple_dir)

    } else {

        ch_purple_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Append RNA data to SAGE VCF
    //
    // channel: [ meta, sage_append_vcf ]
    ch_sage_somatic_append_out = Channel.empty()
    ch_sage_germline_append_out = Channel.empty()
    if (run_config.stages.orange) {

        // NOTE(SW): currently used only for ORANGE but will also be used for Neo once implemented

        SAGE_APPEND(
            ch_inputs,
            ch_purple_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
        )

        ch_versions = ch_versions.mix(SAGE_APPEND.out.versions)
        ch_sage_somatic_append_out = ch_sage_somatic_append_out.mix(SAGE_APPEND.out.somatic_vcf)
        ch_sage_germline_append_out = ch_sage_germline_append_out.mix(SAGE_APPEND.out.germline_vcf)

    } else {

        ch_sage_somatic_append_out = ch_inputs.map { meta -> [meta, []] }
        ch_sage_germline_append_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Group structural variants into higher order events with LINX
    //
    // channel: [ meta, linx_annotation_dir ]
    ch_linx_somatic_out = Channel.empty()
    ch_linx_germline_out = Channel.empty()
    if (run_config.stages.linx) {

        LINX_ANNOTATION(
            ch_inputs,
            ch_purple_out,
            ref_data.genome_version,
            hmf_data.ensembl_data_resources,
            hmf_data.known_fusion_data,
            hmf_data.driver_gene_panel,
            linx_gene_id_file,
        )

        ch_versions = ch_versions.mix(LINX_ANNOTATION.out.versions)

        ch_linx_somatic_out = ch_linx_somatic_out.mix(LINX_ANNOTATION.out.somatic)
        ch_linx_germline_out = ch_linx_germline_out.mix(LINX_ANNOTATION.out.germline)

    } else {

        ch_linx_somatic_out = ch_inputs.map { meta -> [meta, []] }
        ch_linx_germline_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Visualise LINX annotations
    //
    // channel: [ meta, linx_visualiser_dir_reportable ]
    ch_linx_somatic_plot_reportable_out = Channel.empty()
    if (run_config.stages.linx) {

        LINX_PLOTTING(
            ch_inputs,
            ch_linx_somatic_out,
            ref_data.genome_version,
            hmf_data.ensembl_data_resources,
        )

        ch_versions = ch_versions.mix(LINX_PLOTTING.out.versions)

        ch_linx_somatic_plot_reportable_out = ch_linx_somatic_plot_reportable_out.mix(LINX_PLOTTING.out.visualiser_dir_reportable)

    } else {

        ch_linx_somatic_plot_reportable_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run SAMtools flagstat to generate stats required for ORANGE
    //
    // channel: [ meta, metrics ]
    ch_flagstat_somatic_out = Channel.empty()
    ch_flagstat_germline_out = Channel.empty()
    if (run_config.stages.orange && run_config.stages.flagstat) {

        FLAGSTAT_METRICS(
            ch_inputs,
        )

        ch_versions = ch_versions.mix(FLAGSTAT_METRICS.out.versions)

        ch_flagstat_somatic_out = ch_flagstat_somatic_out.mix(FLAGSTAT_METRICS.out.somatic)
        ch_flagstat_germline_out = ch_flagstat_germline_out.mix(FLAGSTAT_METRICS.out.germline)

    } else {

        ch_flagstat_somatic_out = ch_inputs.map { meta -> [meta, []] }
        ch_flagstat_germline_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run Bam Tools to generate stats required for downstream processes
    //
    // channel: [ meta, metrics ]
    ch_bamtools_somatic_out = Channel.empty()
    ch_bamtools_germline_out = Channel.empty()
    if (run_config.stages.bamtools) {

        BAMTOOLS_METRICS(
            ch_inputs,
            ref_data.genome_fasta,
            ref_data.genome_version,
        )

        ch_versions = ch_versions.mix(BAMTOOLS_METRICS.out.versions)

        ch_bamtools_somatic_out = ch_bamtools_somatic_out.mix(BAMTOOLS_METRICS.out.somatic)
        ch_bamtools_germline_out = ch_bamtools_germline_out.mix(BAMTOOLS_METRICS.out.germline)

    } else {

        ch_bamtools_somatic_out = ch_inputs.map { meta -> [meta, []] }
        ch_bamtools_germline_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run Sigs to fit somatic smlv to signature definitions
    //
    // channel: [ meta, sigs_dir ]
    ch_sigs_out = Channel.empty()
    if (run_config.stages.sigs) {

        SIGS_FITTING(
            ch_inputs,
            ch_purple_out,
            hmf_data.sigs_signatures,
        )

        ch_versions = ch_versions.mix(SIGS_FITTING.out.versions)

        ch_sigs_out = ch_sigs_out.mix(SIGS_FITTING.out.sigs_dir)

    } else {

        ch_sigs_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run CHORD to predict HR deficiency status
    //
    // channel: [ meta, chord_dir ]
    ch_chord_out = Channel.empty()
    if (run_config.stages.chord) {

        CHORD_PREDICTION(
            ch_inputs,
            ch_purple_out,
            ref_data.genome_version,
        )

        ch_versions = ch_versions.mix(CHORD_PREDICTION.out.versions)

        ch_chord_out = ch_chord_out.mix(CHORD_PREDICTION.out.chord_dir)

    } else {

        ch_chord_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run LILAC for HLA typing and somatic CNV and SNV calling
    //
    // channel: [ meta, lilac_dir ]
    ch_lilac_out = Channel.empty()
    if (run_config.stages.lilac) {

        // Use HLA slice BED if provided in params or set as default requirement
        ref_data_hla_slice_bed = params.containsKey('ref_data_hla_slice_bed') ? params.ref_data_hla_slice_bed : []

        LILAC_CALLING(
            ch_inputs,
            ch_purple_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            hmf_data.lilac_resources,
            ref_data_hla_slice_bed,
        )

        ch_versions = ch_versions.mix(LILAC_CALLING.out.versions)

        ch_lilac_out = ch_lilac_out.mix(LILAC_CALLING.out.lilac_dir)

    } else {

        ch_lilac_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run VIRUSBreakend and Virus Interpreter to quantify viral content
    //
    // channel: [ meta, virusinterpreter_dir ]
    ch_virusinterpreter_out = Channel.empty()
    if (run_config.stages.virusinterpreter) {

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
            gridss_config,
        )

        ch_versions = ch_versions.mix(VIRUSBREAKEND_CALLING.out.versions)

        ch_virusinterpreter_out = ch_virusinterpreter_out.mix(VIRUSBREAKEND_CALLING.out.virusinterpreter_dir)

    } else {

        ch_virusinterpreter_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run CUPPA predict tissue of origin
    //
    // channel: [ meta, cuppa_dir ]
    ch_cuppa_out = Channel.empty()
    if (run_config.stages.cuppa) {

        CUPPA_PREDICTION(
            ch_inputs,
            ch_isofox_out,
            ch_purple_out,
            ch_linx_somatic_out,
            ch_virusinterpreter_out,
            ref_data.genome_version,
            hmf_data.cuppa_resources,
        )

        ch_versions = ch_versions.mix(CUPPA_PREDICTION.out.versions)

        ch_cuppa_out = ch_cuppa_out.mix(CUPPA_PREDICTION.out.cuppa_dir)

    } else {

        ch_cuppa_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run ORANGE to generate static PDF report
    //
    if (run_config.stages.orange) {

        ORANGE_REPORTING(
            ch_inputs,
            ch_bamtools_somatic_out,
            ch_bamtools_germline_out,
            ch_flagstat_somatic_out,
            ch_flagstat_germline_out,
            ch_sage_somatic_dir_out,
            ch_sage_germline_dir_out,
            ch_sage_somatic_append_out,
            ch_sage_germline_append_out,
            ch_purple_out,
            ch_linx_somatic_out,
            ch_linx_somatic_plot_reportable_out,
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
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
