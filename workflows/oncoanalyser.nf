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
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

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
include { CHECK_SAMPLESHEET } from '../modules/local/check_samplesheet/main'

include { AMBER             } from '../modules/local/amber/main'
include { CHORD             } from '../modules/local/chord/main'
include { COBALT            } from '../modules/local/cobalt/main'
include { CUPPA_CLASSIFIER  } from '../modules/local/cuppa/classifier/main'
include { CUPPA_VISUALISER  } from '../modules/local/cuppa/visualiser/main'
include { ISOFOX            } from '../modules/local/isofox/main'
include { LINX_REPORT       } from '../modules/local/gpgr/linx_report/main'
include { PURPLE            } from '../modules/local/purple/main'
include { SIGS              } from '../modules/local/sigs/main'
include { TEAL              } from '../modules/local/teal/main'
include { VIRUSBREAKEND     } from '../modules/local/virusbreakend/main'
include { VIRUSINTERPRETER  } from '../modules/local/virusinterpreter/main'

//
// SUBWORKFLOWS
//
include { GRIDSS            } from '../subworkflows/local/gridss'
include { GRIPSS            } from '../subworkflows/local/gripss'
include { LILAC             } from '../subworkflows/local/lilac'
include { LINX              } from '../subworkflows/local/linx'
include { PAVE              } from '../subworkflows/local/pave'
include { PREPARE_INPUT     } from '../subworkflows/local/prepare_input'
include { PREPARE_REFERENCE } from '../subworkflows/local/prepare_reference'
include { SAGE              } from '../subworkflows/local/sage'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { PICARD_COLLECTWGSMETRICS as COLLECTWGSMETRICS } from '../modules/nf-core/picard/collectwgsmetrics/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Get absolute file paths
samplesheet   = Utils.getFileObject(params.input)
gridss_config = Utils.getFileObject(params.gridss_config)

workflow ONCOANALYSER {
    // Create channel for versions
    // channel: [versions.yml]
    ch_versions = Channel.empty()

    // Get inputs from samplesheet
    PREPARE_INPUT(
        samplesheet,
    )
    ch_inputs = PREPARE_INPUT.out.data

    // Set up reference data and unpack HMF data map for convenience
    PREPARE_REFERENCE(run)
    hmf_data = PREPARE_REFERENCE.out.hmf_data

    // Set up channel with common inputs for several processes
    if (run.amber || run.cobalt || run.pave || run.lilac || run.teal) {
        // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai]
        ch_bams_and_indices = ch_inputs
            .map { meta ->
                def tumor_bam = meta.get([Constants.FileType.BAM_WGS, Constants.DataType.TUMOR])
                def normal_bam = meta.get([Constants.FileType.BAM_WGS, Constants.DataType.NORMAL])
                [meta, tumor_bam, normal_bam, "${tumor_bam}.bai", "${normal_bam}.bai"]
            }
    }

    // channel (present): [val(meta)]
    // channel (absent): [val(meta)]
    ch_inputs_wts = ch_inputs
        .branch { meta ->
            def key = [Constants.FileType.BAM_WTS, Constants.DataType.TUMOR]
            present: meta.containsKey(key)
                return meta
            absent: ! meta.containsKey(key)
                return meta
        }

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
                def meta_isofox = [key: meta.id, id: meta.get(['sample_name', Constants.DataType.TUMOR])]
                return [meta_isofox, meta.get([Constants.FileType.BAM_WTS, Constants.DataType.TUMOR])]
            }

        // Run process
        ISOFOX(
            ch_isofox_inputs,
            PREPARE_REFERENCE.out.genome_fasta,
            PREPARE_REFERENCE.out.genome_fai,
            PREPARE_REFERENCE.out.genome_version,
            hmf_data.ensembl_data_resources,
            hmf_data.isofox_counts,
            hmf_data.isofox_gc_ratios,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(ISOFOX.out.versions)
        ch_isofox_out = ch_isofox_out.mix(WorkflowOncoanalyser.restoreMeta(ISOFOX.out.isofox_dir, ch_inputs))
    }

    //
    // MODULE: Run COLLECTWGSMETRICS to generate stats required for downstream processes
    //
    if (run.virusinterpreter || run.teal) {
        // Create inputs and create process-specific meta
        // NOTE(SW): CUPPA only requires collectwgsmetrics for the tumor sample in the upstream
        // process Virus Interpreter but TEAL currently requires collectwgsmetrics for both tumor
        // and normal sample
        // channel: [val(meta_cwm), bam]
        ch_cwm_inputs_all = ch_inputs
            .flatMap { meta ->
                def sample_types = run.teal ? [Constants.DataType.TUMOR, Constants.DataType.NORMAL] : [Constants.DataType.TUMOR]
                return sample_types
                    .collect { sample_type ->
                        def bam = meta.get([Constants.FileType.BAM_WGS, sample_type])
                        def sample_name = meta.get(['sample_name', sample_type])
                        def meta_cwm = [
                            key: meta.id,
                            id: sample_name,
                            // NOTE(SW): must use string representation for caching purposes
                            sample_type: sample_type.name(),
                        ]
                        return [meta_cwm, bam]
                    }
            }

        // Collapse duplicate files e.g. repeated normal BAMs for multiple tumor samples
        // NOTE(SW): no effective blocking by .groupTuple() as we're not dependent
        // on any process
        // channel: [val(meta_cwm), bam]
        ch_cwm_inputs = ch_cwm_inputs_all
            .map { [it[1..-1], it[0]] }
            .groupTuple()
            .map { filepaths, meta_cwm ->
                def (keys, sample_names, sample_types) = meta_cwm
                    .collect {
                        [it.key, it.id, it.sample_type]
                    }
                    .transpose()

                def sample_types_unique = sample_types.unique(false)
                assert sample_types_unique.size() == 1
                def sample_type = sample_types_unique[0]

                def meta_cwm_new = [
                    keys: keys,
                    id: sample_names.join('__'),
                    id_simple: keys.join('__'),
                    sample_type: sample_type,
                ]
                return [meta_cwm_new, *filepaths]
            }

        // Run process
        COLLECTWGSMETRICS(
            ch_cwm_inputs,
            PREPARE_REFERENCE.out.genome_fasta,
        )

        // Set outputs, process outputs and restore original meta
        ch_versions = ch_versions.mix(COLLECTWGSMETRICS.out.versions)

        // Replicate outputs to reverse unique operation
        // channel: [val(meta_cwm_individual), sample_type, metrics]
        ch_cwm_output_individual = COLLECTWGSMETRICS.out.metrics
            .flatMap { meta_cwm_shared, metrics ->
                meta_cwm_shared.keys.collect { key ->
                    return [meta_cwm_shared + [key: key], meta_cwm_shared.sample_type, metrics]
                }
            }

        // Match outputs to original meta and set output
        // channel (tumor): [val(meta), metrics]
        // channel (normal): [val(meta), metrics]
        ch_cwm_output = WorkflowOncoanalyser.restoreMeta(ch_cwm_output_individual, ch_inputs)
            .branch { meta, sample_type_str, metrics ->
                def sample_type = Utils.getEnumFromString(sample_type_str, Constants.DataType)
                tumor: sample_type == Constants.DataType.TUMOR
                    return [meta, metrics]
                normal: sample_type == Constants.DataType.NORMAL
                    return [meta, metrics]
            }
    }

    //
    // MODULE: Run AMBER to obtain b-allele frequencies
    //
    // channel: [val(meta), amber_dir]
    ch_amber_out = Channel.empty()
    if (run.amber) {
        // Create inputs and create process-specific meta
        // channel: [val(meta_amber), tumor_bam, normal_bam, tumor_bai, normal_bai]
        ch_amber_inputs = ch_bams_and_indices
            .map {
                def meta = it[0]
                def fps = it[1..-1]
                def meta_amber = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: meta.get(['sample_name', Constants.DataType.TUMOR]),
                    normal_id: meta.get(['sample_name', Constants.DataType.NORMAL]),
                ]
                return [meta_amber, *fps]
            }

        // Run process
        AMBER(
            ch_amber_inputs,
            PREPARE_REFERENCE.out.genome_version,
            hmf_data.heterozygous_sites,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(AMBER.out.versions)
        ch_amber_out = ch_amber_out.mix(WorkflowOncoanalyser.restoreMeta(AMBER.out.amber_dir, ch_inputs))
    }

    //
    // MODULE: Run COBALT to obtain read ratios
    //
    // channel: [val(meta), cobalt_dir]
    ch_cobalt_out = Channel.empty()
    if (run.cobalt) {
        // Create inputs and create process-specific meta
        // channel: [meta_cobalt, tbam, nbam, tbai, nbai]
        ch_cobalt_inputs = ch_bams_and_indices
            .map {
                def meta = it[0]
                def fps = it[1..-1]
                def meta_cobalt = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: meta.get(['sample_name', Constants.DataType.TUMOR]),
                    normal_id: meta.get(['sample_name', Constants.DataType.NORMAL]),
                ]
                return [meta_cobalt, *fps]
            }

        // Run process
        COBALT(
            ch_cobalt_inputs,
            hmf_data.gc_profile,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(COBALT.out.versions)
        ch_cobalt_out = ch_cobalt_out.mix(WorkflowOncoanalyser.restoreMeta(COBALT.out.cobalt_dir, ch_inputs))
    }

    //
    // SUBWORKFLOW: Call structural variants with GRIDSS
    //
    // channel: [val(meta), gridss_vcf]
    ch_gridss_out = Channel.empty()
    if (run.gridss) {
        ch_gridss_inputs = ch_inputs
            .map { meta ->
                def tumor_bam = meta.get([Constants.FileType.BAM_WGS, Constants.DataType.TUMOR])
                def normal_bam = meta.get([Constants.FileType.BAM_WGS, Constants.DataType.NORMAL])
                [meta, tumor_bam, normal_bam]
            }
        GRIDSS(
            ch_gridss_inputs,
            gridss_config,
            PREPARE_REFERENCE.out.genome_fasta,
            PREPARE_REFERENCE.out.genome_fai,
            PREPARE_REFERENCE.out.genome_dict,
            PREPARE_REFERENCE.out.genome_bwa_index,
            PREPARE_REFERENCE.out.genome_bwa_index_image,
            PREPARE_REFERENCE.out.genome_gridss_index,
            hmf_data.gridss_region_blocklist,
        )
        ch_versions = ch_versions.mix(GRIDSS.out.versions)
        ch_gridss_out = ch_gridss_out.mix(GRIDSS.out.results)
    }

    //
    // MODULE: Run GRIPSS to filter GRIDSS SV calls
    //
    // channel: [val(meta), hard_vcf, hard_tbi, soft_vcf, soft_tbi]
    ch_gripss_germline_out = Channel.empty()
    // channel: [val(meta), hard_vcf, hard_tbi, soft_vcf, soft_tbi]
    ch_gripss_somatic_out = Channel.empty()
    if (run.gripss) {
        // Select input source
        // channel: [val(meta), gridss_vcf]
        if (run.gridss) {
            ch_gripss_inputs_source = ch_gridss_out
        } else {
            ch_gripss_inputs_source = WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.GRIDSS_VCF, Constants.DataType.TUMOR_NORMAL])
        }

        // Create inputs and create process-specific meta
        // channel: [val(meta_gripss), gridss_vcf]
        ch_gripss_inputs = ch_gripss_inputs_source
            .map { meta, gridss_vcf ->
                def meta_gripss = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: meta.get(['sample_name', Constants.DataType.TUMOR]),
                    normal_id: meta.get(['sample_name', Constants.DataType.NORMAL]),
                ]
                return [meta_gripss, gridss_vcf]
            }

        // Call subworkflow to run processes
        GRIPSS(
            ch_gripss_inputs,
            PREPARE_REFERENCE.out.genome_fasta,
            PREPARE_REFERENCE.out.genome_fai,
            PREPARE_REFERENCE.out.genome_version,
            hmf_data.gridss_pon_breakends,
            hmf_data.gridss_pon_breakpoints,
            hmf_data.known_fusions,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(GRIPSS.out.versions)
        ch_gripss_germline_out = ch_gripss_germline_out.mix(WorkflowOncoanalyser.restoreMeta(GRIPSS.out.germline, ch_inputs))
        ch_gripss_somatic_out = ch_gripss_somatic_out.mix(WorkflowOncoanalyser.restoreMeta(GRIPSS.out.somatic, ch_inputs))
    }

    //
    // SUBWORKFLOW: call SNV, MNV, and small INDELS with SAGE
    //
    // channel: [val(meta), sage_vcf]
    ch_sage_germline_out = Channel.empty()
    // channel: [val(meta), sage_vcf]
    ch_sage_somatic_out = Channel.empty()
    if (run.sage) {
        // Create inputs and create process-specific meta
        ch_sage_inputs = ch_bams_and_indices
            .map { meta, tbam, nbam, tbai, nbai ->
                def meta_sage = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: meta.get(['sample_name', Constants.DataType.TUMOR]),
                    normal_id: meta.get(['sample_name', Constants.DataType.NORMAL]),
                ]
                return [meta_sage, tbam, nbam, tbai, nbai]
            }

        // Call subworkflow to run processes
        SAGE(
            ch_sage_inputs,
            PREPARE_REFERENCE.out.genome_fasta,
            PREPARE_REFERENCE.out.genome_fai,
            PREPARE_REFERENCE.out.genome_dict,
            PREPARE_REFERENCE.out.genome_version,
            hmf_data.sage_known_hotspots_germline,
            hmf_data.sage_known_hotspots_somatic,
            hmf_data.sage_coding_panel,
            hmf_data.sage_highconf_regions,
            hmf_data.sage_pon,
            hmf_data.segment_mappability,
            hmf_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(SAGE.out.versions)
        ch_sage_germline_out = ch_sage_germline_out.mix(WorkflowOncoanalyser.restoreMeta(SAGE.out.germline, ch_inputs))
        ch_sage_somatic_out = ch_sage_somatic_out.mix(WorkflowOncoanalyser.restoreMeta(SAGE.out.somatic, ch_inputs))
    }

    //
    // SUBWORKFLOW: Annotate variants with PAVE
    //
    // channel: [val(meta), pave_vcf]
    ch_pave_germline_out = Channel.empty()
    // channel: [val(meta), pave_vcf]
    ch_pave_somatic_out = Channel.empty()
    if (run.pave) {
        // Select input sources
        // channel: [meta, sage_vcf]
        if (run.sage) {
            ch_pave_germline_inputs_source = ch_sage_germline_out
            ch_pave_somatic_inputs_source = ch_sage_somatic_out
        } else {
            ch_pave_germline_inputs_source = WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.SAGE_VCF, Constants.DataType.NORMAL])
            ch_pave_somatic_inputs_source = WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.SAGE_VCF, Constants.DataType.TUMOR])
        }

        // Create inputs and create process-specific meta
        // channel: [val(meta_pave), sage_vcf]
        ch_pave_germline_inputs = ch_pave_germline_inputs_source
            .map { meta, sage_vcf ->
                def pave_meta = [
                    key: meta.id,
                    id: meta.get(['sample_name', Constants.DataType.TUMOR]),
                ]
                return [pave_meta, sage_vcf]
            }
        ch_pave_somatic_inputs = ch_pave_somatic_inputs_source
            .map { meta, sage_vcf ->
                def pave_meta = [
                    key: meta.id,
                    id: meta.get(['sample_name', Constants.DataType.TUMOR]),
                ]
                return [pave_meta, sage_vcf]
            }

        // Call subworkflow to run processes
        PAVE(
            ch_pave_germline_inputs,
            ch_pave_somatic_inputs,
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
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(PAVE.out.versions)
        ch_pave_germline_out = ch_pave_germline_out.mix(WorkflowOncoanalyser.restoreMeta(PAVE.out.germline, ch_inputs))
        ch_pave_somatic_out = ch_pave_somatic_out.mix(WorkflowOncoanalyser.restoreMeta(PAVE.out.somatic, ch_inputs))
    }

    //
    // MODULE: Run PURPLE for CNV calling, purity and ploidy inference, SV recovery
    //
    // channel: [val(meta), purple_dir]
    ch_purple_out = Channel.empty()
    if (run.purple) {
        // Select input sources
        // channel: [val(meta), sv_hard_vcf, sv_hard_tbi, sv_soft_vcf, sv_soft_tbi]
        if (run.gripss) {
            ch_purple_inputs_sv = ch_gripss_somatic_out
        } else {
            ch_purple_inputs_sv = ch_inputs
                .map { meta ->
                    def sv_hard_vcf = meta[[Constants.FileType.GRIPSS_HARD_VCF, Constants.DataType.TUMOR]]
                    def sv_soft_vcf = meta[[Constants.FileType.GRIPSS_SOFT_VCF, Constants.DataType.TUMOR]]
                    return [
                        meta,
                        sv_hard_vcf,
                        "${sv_hard_vcf}.tbi",
                        sv_soft_vcf,
                        "${sv_soft_vcf}.tbi",
                    ]
                }
        }

        // channel: [val(meta), amber_dir, cobalt_dir, sv_hard_vcf, sv_hard_tbi, sv_soft_vcf, sv_soft_tbi, smlv_tumor_vcf, smlv_normal_vcf]
        ch_purple_inputs_sources = WorkflowOncoanalyser.groupByMeta(
            run.amber ? ch_amber_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.AMBER_DIR, Constants.DataType.TUMOR_NORMAL]),
            run.cobalt ? ch_cobalt_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.COBALT_DIR, Constants.DataType.TUMOR_NORMAL]),
            ch_purple_inputs_sv,
            run.pave ? ch_pave_somatic_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.PAVE_VCF, Constants.DataType.TUMOR]),
            run.pave ? ch_pave_germline_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.PAVE_VCF, Constants.DataType.NORMAL]),
        )

        // channel: [val(meta_purple), amber_dir, cobalt_dir, sv_hard_vcf, sv_hard_tbi, sv_soft_vcf, sv_soft_tbi, smlv_tumor_vcf, smlv_normal_vcf]
        ch_purple_inputs = ch_purple_inputs_sources
            .map {
                def meta = it[0]
                def other = it[1..-1]
                def meta_purple = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: meta.get(['sample_name', Constants.DataType.TUMOR]),
                    normal_id: meta.get(['sample_name', Constants.DataType.NORMAL]),
              ]
              return [meta_purple, *other]
            }

        PURPLE(
            ch_purple_inputs,
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
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(PURPLE.out.versions)
        ch_purple_out = ch_purple_out.mix(WorkflowOncoanalyser.restoreMeta(PURPLE.out.purple_dir, ch_inputs))
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
            ch_sigs_inputs_source = WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.PURPLE_DIR, Constants.DataType.TUMOR_NORMAL])
        }

        // Create inputs and create process-specific meta
        // channel: [val(meta_sigs), smlv_vcf]
        ch_sigs_inputs = ch_sigs_inputs_source
            .map { meta, purple_dir ->
                def smlv_vcf = file(purple_dir).resolve("${meta.get(['sample_name', Constants.DataType.TUMOR])}.purple.somatic.vcf.gz")
                def meta_sigs = [
                    id: meta.id,
                    tumor_id: meta.get(['sample_name', Constants.DataType.TUMOR]),
                ]
                return [meta_sigs, smlv_vcf]
            }

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
    if (run.chord) {
        // Select input sources
        // channel: [val(meta), purple_dir]
        if (run.purple) {
          ch_chord_inputs_source = ch_purple_out
        } else {
          ch_chord_inputs_source = WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.PURPLE_DIR, Constants.DataType.TUMOR_NORMAL])
        }

        // Create inputs and create process-specific meta
        // channel: [val(meta), smlv_vcf, sv_vcf]
        ch_chord_inputs = ch_chord_inputs_source
            .map { meta, purple_dir ->
                def tumor_id = meta.get(['sample_name', Constants.DataType.TUMOR])
                def smlv_vcf = file(purple_dir).resolve("${tumor_id}.purple.somatic.vcf.gz")
                def sv_vcf = file(purple_dir).resolve("${tumor_id}.purple.sv.vcf.gz")
                def meta_chord = [id: meta.id]
                return [meta_chord, smlv_vcf, sv_vcf]
            }

        // Run process
        CHORD(
          ch_chord_inputs,
          PREPARE_REFERENCE.out.genome_version,
        )

        // Set outputs
        ch_versions = ch_versions.mix(CHORD.out.versions)
    }

    //
    // MODULE: Run TEAL to characterise teleomeres
    //
    if (run.teal) {
        // Select input sources
        // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai, tumor_wgs_metrics, normal_wgs_metrics, cobalt_dir, purple_dir]
        // NOTE(SW): assuming here that TEAL is being run in tumor/normal mode and so we expect a tumor metrics file and normal metrics file
        ch_teal_inputs_source = WorkflowOncoanalyser.groupByMeta(
            ch_bams_and_indices,
            run.collectwgsmetrics ? ch_cwm_output.tumor : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.COLLECTWGSMETRICS, Constants.DataType.TUMOR]),
            run.collectwgsmetrics ? ch_cwm_output.normal : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.COLLECTWGSMETRICS, Constants.DataType.NORMAL]),
            run.cobalt ? ch_cobalt_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.COBALT_DIR, Constants.DataType.TUMOR_NORMAL]),
            run.purple ? ch_purple_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.PURPLE_DIR, Constants.DataType.TUMOR_NORMAL]),
        )

        // Create inputs and create process-specific meta
        // channel: [val(meta_teal), tumor_bam, normal_bam, tumor_bai, normal_bai, tumor_wgs_metrics, normal_wgs_metrics, cobalt_dir, purple_dir]
        ch_teal_inputs = ch_teal_inputs_source
            .map {
                def meta = it[0]
                def other = it[1..-1]
                def meta_teal = [
                    id: meta.id,
                    tumor_id: meta.get(['sample_name', Constants.DataType.TUMOR]),
                    normal_id: meta.get(['sample_name', Constants.DataType.NORMAL]),
                ]
                return [meta_teal, *other]
            }

        // Run process
        TEAL(
            ch_teal_inputs,
        )

        // Set outputs
        ch_versions = ch_versions.mix(TEAL.out.versions)
    }

    //
    // SUBWORKFLOW: Run LILAC for HLA typing and somatic CNV and SNV calling
    //
    if (run.lilac) {
        // Select input sources
        // channel: [val(meta), purple_dir]
        if (run.purple) {
            ch_lilac_inputs_source = ch_purple_out
        } else {
            ch_lilac_inputs_source = WorkflowOncoanalyser.getInput([Constants.FileType.PURPLE_DIR, Constants.DataType.TUMOR_NORMAL])
        }

        // Call subworkflow to run processes
        LILAC(
            ch_bams_and_indices,
            run.purple ? ch_purple_out : WorkflowOncoanalyser.getInput([Constants.FileType.PURPLE_DIR, Constants.DataType.TUMOR_NORMAL]),
            PREPARE_REFERENCE.out.genome_fasta,
            PREPARE_REFERENCE.out.genome_fai,
            hmf_data.lilac_resources,
        )

        // Set outputs
        ch_versions = ch_versions.mix(LILAC.out.versions)
    }

    //
    // MODULE: Run VIRUSBreakend and Virus Interpreter to quantify viral content
    //
    // NOTE(SW): kept separate from CUPPA conditional block since we'll allow users to run this independently
    ch_virusinterpreter_out = Channel.empty()
    if (run.virusinterpreter) {
        // Create inputs and create process-specific meta
        // channel: [val(meta_virus), tumor_bam]
        ch_virusbreakend_inputs = ch_inputs
            .map { meta ->
                def meta_virus = [
                    key: meta.id,
                    id: meta.id,
                ]
                return [meta_virus, meta.get([Constants.FileType.BAM_WGS, Constants.DataType.TUMOR])]
            }

        // Run process
        VIRUSBREAKEND(
            ch_virusbreakend_inputs,
            gridss_config,
            PREPARE_REFERENCE.out.genome_fasta,
            PREPARE_REFERENCE.out.genome_fai,
            PREPARE_REFERENCE.out.genome_dict,
            PREPARE_REFERENCE.out.genome_bwa_index,
            PREPARE_REFERENCE.out.genome_bwa_index_image,
            PREPARE_REFERENCE.out.genome_gridss_index,
            PREPARE_REFERENCE.out.virusbreakenddb,
        )

        // Set outputs
        ch_versions = ch_versions.mix(VIRUSBREAKEND.out.versions)

        // Create inputs and create process-specific meta
        // channel: [val(meta), purple_purity, purple_qc]
        ch_virusinterpreter_inputs_purple = ch_purple_out
            .map { meta, purple_dir ->
                def purple_purity = file(purple_dir).resolve("${meta.get(['sample_name', Constants.DataType.TUMOR])}.purple.purity.tsv")
                def purple_qc = file(purple_dir).resolve("${meta.get(['sample_name', Constants.DataType.TUMOR])}.purple.qc")
                return [meta, purple_purity, purple_qc]
            }

        // channel: [val(meta), virus_tsv, purple_purity, purple_qc, wgs_metrics]
        ch_virusinterpreter_inputs_full = WorkflowOncoanalyser.groupByMeta(
            WorkflowOncoanalyser.restoreMeta(VIRUSBREAKEND.out.tsv, ch_inputs),
            ch_virusinterpreter_inputs_purple,
            ch_cwm_output.tumor,
        )

        // Create inputs and create process-specific meta
        // channel: [val(meta_virus), virus_tsv, purple_purity, purple_qc, wgs_metrics]
        ch_virusinterpreter_inputs = ch_virusinterpreter_inputs_full
            .map {
                def meta = it[0]
                def meta_virus = [
                    key: meta.id,
                    id: meta.get(['sample_name', Constants.DataType.TUMOR]),
                ]
                return [meta_virus, *it[1..-1]]
            }

        // Run process
        VIRUSINTERPRETER(
            ch_virusinterpreter_inputs,
            hmf_data.virus_taxonomy_db,
            hmf_data.virus_reporting_db,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(VIRUSINTERPRETER.out.versions)
        ch_virusinterpreter_out = ch_virusinterpreter_out.mix(WorkflowOncoanalyser.restoreMeta(VIRUSINTERPRETER.out.virusinterpreter_dir, ch_inputs))
    }

    //
    // SUBWORKFLOW: Group structural variants into higher order events with LINX
    //
    // channel: [val(meta), linx_annotation_dir, linx_visuliaser_dir]
    ch_linx_somatic_out = Channel.empty()
    if (run.linx) {
        // Select input sources
        // channel: [val(meta), vcf]
        if (run.gripss) {
            ch_linx_inputs_germline_source = ch_gripss_germline_out.map { meta, h, hi, s, si -> [meta, h] }
        } else {
            ch_linx_inputs_germline_source = WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.GRIPSS_HARD_VCF, Constants.DataType.NORMAL])
        }
        // channel: [val(meta), sv_vcf, purple_dir]
        ch_linx_inputs_somatic_source = run.purple ? ch_purple_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.PURPLE_DIR, Constants.DataType.TUMOR_NORMAL])

        // Create inputs and create process-specific meta
        // channel: [val(meta_linx), sv_vcf]
        ch_linx_inputs_germline = ch_linx_inputs_germline_source
            .map {
                def meta = it[0]
                def meta_linx = [
                    key: meta.id,
                    id: meta.get(['sample_name', Constants.DataType.NORMAL]),
                ]
                return [meta_linx, it[1..-1]]
            }
        // channel: [val(meta_linx), purple_dir]
        ch_linx_inputs_somatic = ch_linx_inputs_somatic_source
            .map {
                def meta = it[0]
                def meta_linx = [
                    key: meta.id,
                    id: meta.get(['sample_name', Constants.DataType.TUMOR]),
                ]
                return [meta_linx, it[1..-1]]
            }

        // Run process
        LINX(
            ch_linx_inputs_germline,
            ch_linx_inputs_somatic,
            PREPARE_REFERENCE.out.genome_version,
            hmf_data.linx_fragile_regions,
            hmf_data.linx_lines,
            hmf_data.ensembl_data_resources,
            hmf_data.known_fusion_data,
            hmf_data.driver_gene_panel,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(LINX.out.versions)
        ch_linx_somatic_out = ch_linx_somatic_out.mix(WorkflowOncoanalyser.restoreMeta(LINX.out.somatic, ch_inputs))

        //
        // MODULE: Run gpgr to generate a LINX report
        //
        // Create inputs and create process-specific meta
        // channel: [meta(meta_linx_report), anno_dir, vis_dir]
        ch_linx_report_inputs = ch_linx_somatic_out
            .map { meta, anno_dir, vis_dir ->
                def meta_linx_report = [
                    key: meta.id,
                    id: meta.get(['sample_name', Constants.DataType.TUMOR]),
                ]
                return [meta_linx_report, anno_dir, vis_dir]
            }

        // Run process
        LINX_REPORT(
            ch_linx_report_inputs,
        )

        // Set outputs
        ch_versions = ch_versions.mix(LINX_REPORT.out.versions)
    }

    //
    // MODULE: Run CUPPA predict tissue of origin
    //
    if (run.cuppa) {
        // Select input sources
        // channel: [val(meta), isofox_dir]
        if (run.isofox) {
            ch_cuppa_inputs_isofox = Channel.empty()
                .mix(
                    ch_inputs_wts.absent.map { meta -> [meta, []] },
                    ch_isofox_out,
                )
        } else {
            ch_cuppa_inputs_isofox = ch_inputs
                .map { meta -> [meta, meta.get([Constants.FileType.ISOFOX_DIR, Constants.DataType.TUMOR], [])] }
        }

        // channel: [val(meta), linx_somatic_annotation_dir]
        ch_linx_anno = ch_linx_somatic_out.map { meta, anno_dir, vis_dir -> [meta, anno_dir]}

        // channel: [val(meta), isofox_dir, purple_dir, linx_dir, virusinterpreter_dir]
        // NOTE(SW): the Groovy Collection.flatten method used in
        // WorkflowOncoanalyser.groupByMeta removes optional Isofox input; flattening
        // done manually below to preserve
        ch_cuppa_inputs_source = WorkflowOncoanalyser.groupByMeta(
            ch_cuppa_inputs_isofox,
            run.purple ? ch_purple_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.PURPLE_DIR, Constants.DataType.TUMOR_NORMAL]),
            run.linx ? ch_linx_anno : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.LINX_DIR, Constants.DataType.TUMOR_NORMAL]),
            run.virusinterpreter ? ch_virusinterpreter_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.VIRUSINTERPRETER, Constants.DataType.TUMOR]),
            flatten: false,
        )
            .map { data ->
                def meta = data[0]
                def inputs = data[1..-1].collect { it[0] }
                return [meta, *inputs]
            }

        // Create inputs and create process-specific meta
        // channel: [val(meta_cuppa), isofox_dir, purple_dir, linx_dir, virusinterpreter_dir]
        ch_cuppa_inputs = ch_cuppa_inputs_source
            .map {
                def meta = it[0]
                def meta_cuppa = [
                    key: meta.id,
                    id: meta.get(['sample_name', Constants.DataType.TUMOR]),
                ]
                return [meta_cuppa, *it[1..-1]]
            }

        // Run process
        CUPPA_CLASSIFIER(
            ch_cuppa_inputs,
            hmf_data.cuppa_resources,
        )

        // Set outputs
        ch_versions = ch_versions.mix(CUPPA_CLASSIFIER.out.versions)

        // Run process
        CUPPA_VISUALISER(
            CUPPA_CLASSIFIER.out.csv,
        )

        // Set outputs
        ch_versions = ch_versions.mix(CUPPA_VISUALISER.out.versions)
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
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
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
