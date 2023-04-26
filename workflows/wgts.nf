import Constants
import Processes
import Utils


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowOncoanalyser.initialise(params, workflow, log)

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
// MODULES
//
include { CHORD            } from '../modules/local/chord/main'
include { ISOFOX           } from '../modules/local/isofox/main'
include { PEACH            } from '../modules/local/peach/main'
include { PROTECT          } from '../modules/local/protect/main'
include { SIGS             } from '../modules/local/sigs/main'

//
// SUBWORKFLOWS
//
include { AMBER_PROFILING       } from '../subworkflows/local/amber_profiling'
include { BAMTOOLS_METRICS      } from '../subworkflows/local/bamtools_metrics'
include { COBALT_PROFILING      } from '../subworkflows/local/cobalt_profiling'
include { CUPPA_PREDICTION      } from '../subworkflows/local/cuppa_prediction'
include { GRIDSS_CALLING        } from '../subworkflows/local/gridss_calling'
include { GRIDSS_SVPREP_CALLING } from '../subworkflows/local/gridss_svprep_calling'
include { GRIPSS_FILTERING      } from '../subworkflows/local/gripss_filtering'
include { LILAC_CALLING         } from '../subworkflows/local/lilac_calling'
include { LINX_ANNOTATION       } from '../subworkflows/local/linx_annotation'
include { LINX_PLOTTING         } from '../subworkflows/local/linx_plotting'
include { ORANGE_REPORTING      } from '../subworkflows/local/orange_reporting'
include { PAVE_ANNOTATION       } from '../subworkflows/local/pave_annotation'
include { PREPARE_INPUT         } from '../subworkflows/local/prepare_input'
include { PREPARE_REFERENCE     } from '../subworkflows/local/prepare_reference'
include { PURPLE_CALLING        } from '../subworkflows/local/purple_calling'
include { SAGE_CALLING          } from '../subworkflows/local/sage_calling'
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

    // Get inputs from samplesheet
    // channel: [val(meta)]
    PREPARE_INPUT(
        samplesheet,
    )
    ch_inputs = PREPARE_INPUT.out.data

    // Split inputs WTS and WGS
    // NOTE(SW): assuming there are only t/n pairs i.e. no tumor-only or normal-only

    // channel (present): [val(meta)]
    // channel (absent): [val(meta)]
    ch_inputs_wgs = ch_inputs
        .branch { meta ->
            def key_tumor = [Constants.FileType.BAM, Constants.SampleType.TUMOR, Constants.SequenceType.WGS]
            def key_normal = [Constants.FileType.BAM, Constants.SampleType.NORMAL, Constants.SequenceType.WGS]
            present: meta.containsKey(key_tumor) && meta.containsKey(key_normal)
                return meta
            absent: true
                return meta
        }

    // channel (present): [val(meta)]
    // channel (absent): [val(meta)]
    ch_inputs_wts = ch_inputs
        .branch { meta ->
            def key = [Constants.FileType.BAM, Constants.SampleType.TUMOR, Constants.SequenceType.WTS]
            present: meta.containsKey(key)
                return meta
            absent: ! meta.containsKey(key)
                return meta
        }

    // Set up reference data and unpack HMF data map for convenience
    PREPARE_REFERENCE(run)
    hmf_data = PREPARE_REFERENCE.out.hmf_data

    //
    // MODULE: Run Isofox to analyse WTS data
    //
    // channel: [meta, isofox_dir]
    ch_isofox_out = Channel.empty()
    if (run.isofox) {
        // Create inputs and create process-specific meta
        // channel: [meta_isofox, tumor_bam_wts]
        ch_isofox_inputs = ch_inputs_wts.present
            .map { meta ->
                def bam = Utils.getTumorWtsBam(meta)
                def meta_isofox = [key: meta.id, id: Utils.getTumorWtsSampleName(meta)]
                return [meta_isofox, bam, "${bam}.bai"]
            }

        // Set Isofox cache files
        // NOTE(SW): the Isofox expected count file is read length dependent so required users to explicitly use expect
        // counts generated for 151 bp reads that is available in the HMF reference bundle. When not specifying an
        // expected count file, Isofox will automatically create one for the computed read length. However, doing so
        // greatly increases runtime.
        // NOTE(SW): consider alternative approaches for using the expected count file e.g. generate once at runtime,
        // then use for all samples; generate all possible read lengths outside of pipeline and store on a remote for
        // retrieval at runtime (requires inference of read length)

        // TODO(SW): this must be improved to allow users to set input file, use cache, or generate at runtime;
        // currently does not update functions
        // NOTE(SW): forcing use of cache for now since this feature is incomplete

        //isofox_counts = params.use_isofox_exp_counts_cache ? hmf_data.isofox_counts : []
        isofox_counts = hmf_data.isofox_counts

        // Run process
        ISOFOX(
            ch_isofox_inputs,
            params.isofox_functions,
            PREPARE_REFERENCE.out.genome_fasta,
            PREPARE_REFERENCE.out.genome_fai,
            PREPARE_REFERENCE.out.genome_version,
            hmf_data.ensembl_data_resources,
            isofox_counts,
            hmf_data.isofox_gc_ratios,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(ISOFOX.out.versions)
        ch_isofox_out = ch_isofox_out.mix(WorkflowOncoanalyser.restoreMeta(ISOFOX.out.isofox_dir, ch_inputs))
    }

    //
    // MODULE: Run Bam Tools to generate stats required for downstream processes
    //
    // channel: [val(meta), metrics]
    ch_bamtools_somatic_out = Channel.empty()
    ch_bamtools_germline_out = Channel.empty()
    if (run.bamtools) {

        BAMTOOLS_METRICS(
            ch_inputs_wgs.present,
            PREPARE_REFERENCE.out.genome_fasta,
            PREPARE_REFERENCE.out.genome_version,
            run,
        )

        ch_versions = ch_versions.mix(BAMTOOLS_METRICS.out.versions)
        ch_bamtools_somatic_out = ch_bamtools_somatic_out.mix(BAMTOOLS_METRICS.out.somatic)
        ch_bamtools_germline_out = ch_bamtools_germline_out.mix(BAMTOOLS_METRICS.out.germline)
    }

    //
    // MODULE: Run AMBER to obtain b-allele frequencies
    //
    // channel: [val(meta), amber_dir]
    ch_amber_out = Channel.empty()
    if (run.amber) {

        AMBER_PROFILING(
            ch_inputs_wgs.present,
            PREPARE_REFERENCE.out.genome_version,
            hmf_data.heterozygous_sites,
        )

        ch_versions = ch_versions.mix(AMBER_PROFILING.out.versions)
        ch_amber_out = ch_amber_out.mix(AMBER_PROFILING.out.amber_dir)
    }

    //
    // MODULE: Run COBALT to obtain read ratios
    //
    // channel: [val(meta), cobalt_dir]
    ch_cobalt_out = Channel.empty()
    if (run.cobalt) {

        COBALT_PROFILING(
            ch_inputs_wgs.present,
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
                ch_inputs_wgs.present,
                gridss_config,
                PREPARE_REFERENCE.out.genome_fasta,
                PREPARE_REFERENCE.out.genome_version,
                PREPARE_REFERENCE.out.genome_fai,
                PREPARE_REFERENCE.out.genome_dict,
                PREPARE_REFERENCE.out.genome_bwa_index,
                PREPARE_REFERENCE.out.genome_bwa_index_image,
                PREPARE_REFERENCE.out.genome_gridss_index,
                hmf_data.gridss_region_blocklist,
                hmf_data.sv_prep_blocklist,
                hmf_data.known_fusions,
            )

            ch_versions = ch_versions.mix(GRIDSS_SVPREP_CALLING.out.versions)
            ch_gridss_out = ch_gridss_out.mix(GRIDSS_SVPREP_CALLING.out.results)

        } else {

            GRIDSS_CALLING(
                ch_inputs_wgs.present,
                gridss_config,
                PREPARE_REFERENCE.out.genome_fasta,
                PREPARE_REFERENCE.out.genome_fai,
                PREPARE_REFERENCE.out.genome_dict,
                PREPARE_REFERENCE.out.genome_bwa_index,
                PREPARE_REFERENCE.out.genome_bwa_index_image,
                PREPARE_REFERENCE.out.genome_gridss_index,
                hmf_data.gridss_region_blocklist,
            )

            ch_versions = ch_versions.mix(GRIDSS_CALLING.out.versions)
            ch_gridss_out = ch_gridss_out.mix(GRIDSS_CALLING.out.results)
        }
    }

    //
    // MODULE: Run GRIPSS to filter GRIDSS SV calls
    //
    // channel: [val(meta), vcf, tbi]
    ch_gripss_somatic_out = Channel.empty()
    ch_gripss_somatic_unfiltered_out = Channel.empty()
    ch_gripss_germline_out = Channel.empty()
    if (run.gripss) {

        GRIPSS_FILTERING(
            ch_inputs,
            ch_gridss_out,
            PREPARE_REFERENCE.out.genome_fasta,
            PREPARE_REFERENCE.out.genome_fai,
            PREPARE_REFERENCE.out.genome_version,
            hmf_data.gridss_pon_breakends,
            hmf_data.gridss_pon_breakpoints,
            hmf_data.known_fusions,
            hmf_data.repeatmasker_annotations,
            run,
        )

        ch_versions = ch_versions.mix(GRIPSS_FILTERING.out.versions)
        ch_gripss_somatic_out = ch_gripss_somatic_out.mix(GRIPSS_FILTERING.out.somatic)
        ch_gripss_somatic_unfiltered_out = ch_gripss_somatic_unfiltered_out.mix(GRIPSS_FILTERING.out.somatic_unfiltered)
        ch_gripss_germline_out = ch_gripss_germline_out.mix(GRIPSS_FILTERING.out.germline)
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
            ch_inputs_wgs.present,
            PREPARE_REFERENCE.out.genome_fasta,
            PREPARE_REFERENCE.out.genome_fai,
            PREPARE_REFERENCE.out.genome_dict,
            PREPARE_REFERENCE.out.genome_version,
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
            PREPARE_REFERENCE.out.genome_fasta,
            PREPARE_REFERENCE.out.genome_fai,
            PREPARE_REFERENCE.out.genome_version,
            hmf_data.sage_pon,
            hmf_data.sage_blocklist_regions,
            hmf_data.sage_blocklist_sites,
            hmf_data.clinvar_annotations,
            hmf_data.segment_mappability,
            hmf_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
            hmf_data.gnomad_pon_dir,
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
            ch_gripss_somatic_unfiltered_out,
            PREPARE_REFERENCE.out.genome_fasta,
            PREPARE_REFERENCE.out.genome_fai,
            PREPARE_REFERENCE.out.genome_dict,
            PREPARE_REFERENCE.out.genome_version,
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
    // MODULE: Run Sigs to fit somatic smlv to signature definitions
    //
    if (run.sigs) {
        // Select input sources
        // channel: [val(meta), purple_dir]
        if (run.purple) {
            ch_sigs_inputs_source = ch_purple_out
        } else {
            ch_sigs_inputs_source = WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PURPLE_DIR)
        }

        // Create inputs and create process-specific meta
        // channel: [val(meta_sigs), smlv_vcf]
        ch_sigs_inputs = ch_sigs_inputs_source
            .map { meta, purple_dir ->
                def tumor_id = Utils.getTumorWgsSampleName(meta)
                def smlv_vcf = file(purple_dir).resolve("${tumor_id}.purple.somatic.vcf.gz")

                // Require smlv VCF from the PURPLE directory
                if (!smlv_vcf.exists()) {
                    return Constants.META_PLACEHOLDER
                }

                def meta_sigs = [
                    id: meta.id,
                    tumor_id: meta.getAt(['sample_name', Constants.SampleType.TUMOR]),
                ]
                return [meta_sigs, smlv_vcf]
            }
            .filter { it[0] != Constants.META_PLACEHOLDER }

        SIGS(
          ch_sigs_inputs,
          hmf_data.sigs_signatures,
        )

        // Set outputs
        ch_versions = ch_versions.mix(SIGS.out.versions)
    }

    //
    // MODULE: Run CHORD to predict HR deficiency status
    //
    // channel: [val(meta), chord_prediction]
    ch_chord_out = Channel.empty()
    if (run.chord) {
        // Select input sources
        // channel: [val(meta), purple_dir]
        if (run.purple) {
          ch_chord_inputs_source = ch_purple_out
        } else {
          ch_chord_inputs_source = WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PURPLE_DIR)
        }

        // Create inputs and create process-specific meta
        // channel: [val(meta), smlv_vcf, sv_vcf]
        ch_chord_inputs = ch_chord_inputs_source
            .map { meta, purple_dir ->
                def tumor_id = Utils.getTumorWgsSampleName(meta)
                def smlv_vcf = file(purple_dir).resolve("${tumor_id}.purple.somatic.vcf.gz")
                def sv_vcf = file(purple_dir).resolve("${tumor_id}.purple.sv.vcf.gz")

                // Require both SV and smlv VCF from the PURPLE directory
                if (!smlv_vcf.exists() || !sv_vcf.exists()) {
                    return Constants.META_PLACEHOLDER
                }

                def meta_chord = [key: meta.id, id: meta.id]
                return [meta_chord, smlv_vcf, sv_vcf]
            }
            .filter { it[0] != Constants.META_PLACEHOLDER }

        // Run process
        CHORD(
          ch_chord_inputs,
          PREPARE_REFERENCE.out.genome_version,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(CHORD.out.versions)
        ch_chord_out = ch_chord_out.mix(WorkflowOncoanalyser.restoreMeta(CHORD.out.prediction, ch_inputs))
    }

    //
    // SUBWORKFLOW: Run LILAC for HLA typing and somatic CNV and SNV calling
    //
    // channel: [val(meta), lilac_dir]
    ch_lilac_out = Channel.empty()
    if (run.lilac) {

        // TODO(SW): improve interface here
        LILAC_CALLING(
            ch_inputs,
            ch_inputs_wgs.present,
            ch_inputs_wgs.absent,
            ch_inputs_wts.present,
            ch_inputs_wts.absent,
            ch_purple_out,
            PREPARE_REFERENCE.out.genome_fasta,
            PREPARE_REFERENCE.out.genome_fai,
            hmf_data.lilac_resources,
            run,
        )

        ch_versions = ch_versions.mix(LILAC_CALLING.out.versions)
        ch_lilac_out = ch_lilac_out.mix(LILAC_CALLING.out.lilac_dir)
    }

    //
    // SUBWORKFLOW: Run VIRUSBreakend and Virus Interpreter to quantify viral content
    //
    ch_virusinterpreter_out = Channel.empty()
    if (run.virusinterpreter) {

        VIRUSBREAKEND_CALLING(
            ch_inputs_wgs.present,
            ch_purple_out,
            ch_bamtools_somatic_out,
            PREPARE_REFERENCE.out.genome_fasta,
            PREPARE_REFERENCE.out.genome_fai,
            PREPARE_REFERENCE.out.genome_dict,
            PREPARE_REFERENCE.out.genome_bwa_index,
            PREPARE_REFERENCE.out.genome_bwa_index_image,
            PREPARE_REFERENCE.out.genome_gridss_index,
            PREPARE_REFERENCE.out.virusbreakenddb,
            hmf_data.virus_taxonomy_db,
            hmf_data.virus_reporting_db,
            run,
            gridss_config,
        )

        ch_versions = ch_versions.mix(VIRUSBREAKEND_CALLING.out.versions)
        ch_virusinterpreter_out = ch_virusinterpreter_out.mix(VIRUSBREAKEND_CALLING.out.virusinterpreter)
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
            ch_gripss_germline_out,
            ch_purple_out,
            PREPARE_REFERENCE.out.genome_version,
            hmf_data.linx_fragile_regions,
            hmf_data.linx_lines,
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
            PREPARE_REFERENCE.out.genome_version,
            hmf_data.ensembl_data_resources,
        )

        ch_versions = ch_versions.mix(LINX_PLOTTING.out.versions)
        ch_linx_somatic_plot_out = ch_linx_somatic_plot_out.mix(LINX_PLOTTING.out.visualiser_dir)
    }

    //
    // MODULE: Run PROTECT to match somatic genomic features with treatment evidence
    //
    // channel: [val(meta), protect]
    ch_protect_out = Channel.empty()
    if (run.protect) {
        // Select input sources
        // channel: [val(meta), chord_prediction, purple_dir, linx_dir, virusinterpreter]
        ch_protect_inputs_source = WorkflowOncoanalyser.groupByMeta(
            run.chord ? ch_chord_out : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.CHORD_PREDICTION),
            run.purple ? ch_purple_out : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PURPLE_DIR),
            run.linx ? ch_linx_somatic_out : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.LINX_ANNO_DIR_TUMOR),
            run.lilac ? ch_lilac_out : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.LILAC_DIR),
            run.virusinterpreter ? ch_virusinterpreter_out : WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.VIRUSINTERPRETER_TSV),
        )

        // Create process-specific meta
        // channel: [val(meta_protect), chord_prediction, purple_dir, linx_dir, virusinterpreter]
        ch_protect_inputs = ch_protect_inputs_source
            .map {
                def meta = it[0]
                def other = it[1..-1]
                def meta_protect = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: Utils.getTumorWgsSampleName(meta),
                    normal_id: Utils.getNormalWgsSampleName(meta),
              ]
              return [meta_protect, *other]
            }

        // Run process
        PROTECT(
          ch_protect_inputs,
          PREPARE_REFERENCE.out.genome_version,
          hmf_data.serve_resources,
          hmf_data.driver_gene_panel,
          hmf_data.disease_ontology,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(PROTECT.out.versions)
        ch_protect_out = ch_protect_out.mix(WorkflowOncoanalyser.restoreMeta(PROTECT.out.tsv, ch_inputs))
    }

    //
    // MODULE: Run PEACH to match germline SNVs with pharmacogenetic evidence
    //
    // channel: [val(meta), peach_genotype]
    ch_peach_out = Channel.empty()
    if (run.peach) {
        // Select input sources
        // channel: [val(meta), purple_dir]
        if (run.purple) {
            ch_peach_inputs_source = ch_purple_out
        } else {
            ch_peach_inputs_source = WorkflowOncoanalyser.getInput(ch_inputs, Constants.INPUT.PURPLE_DIR)
        }

        // Create inputs and create process-specific meta
        // channel: [meta_peach, purple_germline_vcf]
        ch_peach_inputs = ch_peach_inputs_source
            .map { meta, purple_dir ->
                def meta_peach = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: Utils.getTumorWgsSampleName(meta),
                    normal_id: Utils.getNormalWgsSampleName(meta),
                ]
                def tumor_id = Utils.getTumorWgsSampleName(meta)
                def purple_germline_vcf = file(purple_dir).resolve("${tumor_id}.purple.germline.vcf.gz")

                if (!purple_germline_vcf.exists()) {
                    return Constants.META_PLACEHOLDER
                }

                return [meta_peach, purple_germline_vcf]
            }
            .filter { it[0] != Constants.META_PLACEHOLDER }

        // Run process
        PEACH(
            ch_peach_inputs,
            PREPARE_REFERENCE.out.genome_version,
            hmf_data.peach_panel,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(PEACH.out.versions)
        ch_peach_out = ch_peach_out.mix(WorkflowOncoanalyser.restoreMeta(PEACH.out.genotype, ch_inputs))
    }

    //
    // SUBWORKFLOW: Run CUPPA predict tissue of origin
    //
    // channel: [val(meta), cuppa_results]
    ch_cuppa_out = Channel.empty()
    // channel: [val(meta), cuppa_summary_plot]
    ch_cuppa_summary_plot_out = Channel.empty()
    // channel: [val(meta), cuppa_feature_plot]
    ch_cuppa_feature_plot_out = Channel.empty()
    if (run.cuppa) {

        CUPPA_PREDICTION(
            ch_inputs,
            ch_inputs_wts.absent,
            ch_isofox_out,
            ch_purple_out,
            ch_linx_somatic_out,
            ch_virusinterpreter_out,
            PREPARE_REFERENCE.out.genome_version,
            hmf_data.cuppa_resources,
            run,
        )

        ch_versions = ch_versions.mix(CUPPA_PREDICTION.out.versions)
        ch_cuppa_out = ch_cuppa_out.mix(CUPPA_PREDICTION.out.csv)
        ch_cuppa_summary_plot_out = ch_cuppa_summary_plot_out.mix(CUPPA_PREDICTION.out.summary_plot)
        ch_cuppa_feature_plot_out = ch_cuppa_feature_plot_out.mix(CUPPA_PREDICTION.out.feature_plot)
    }

    //
    // SUBWORKFLOW: Run ORANGE to generate static PDF report
    //
    if (run.orange) {

        ORANGE_REPORTING(
            ch_inputs,
            ch_inputs_wgs.present,
            ch_bamtools_somatic_out,
            ch_bamtools_germline_out,
            ch_chord_out,
            ch_lilac_out,
            ch_sage_somatic_tumor_bqr_out,
            ch_sage_somatic_normal_bqr_out,
            ch_sage_germline_coverage_out,
            ch_purple_out,
            ch_linx_somatic_out,
            ch_linx_somatic_plot_out,
            ch_linx_germline_out,
            ch_protect_out,
            ch_peach_out,
            ch_cuppa_out,
            ch_cuppa_feature_plot_out,
            ch_cuppa_summary_plot_out,
            ch_virusinterpreter_out,
            PREPARE_REFERENCE.out.genome_version,
            hmf_data.disease_ontology,
            hmf_data.known_fusion_data,
            hmf_data.driver_gene_panel,
            hmf_data.cohort_mapping,
            hmf_data.cohort_percentiles,
            run,
        )

        ch_versions = ch_versions.mix(ORANGE_REPORTING.out.versions)
    }

    //
    // MODULE: Pipeline reporting
    //
    // Run process
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
