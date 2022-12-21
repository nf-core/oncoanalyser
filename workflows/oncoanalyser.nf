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
include { AMBER             } from '../modules/local/amber/main'
include { CHORD             } from '../modules/local/chord/main'
include { COBALT            } from '../modules/local/cobalt/main'
include { CUPPA_CLASSIFIER  } from '../modules/local/cuppa/classifier/main'
include { CUPPA_VISUALISER  } from '../modules/local/cuppa/visualiser/main'
include { ISOFOX            } from '../modules/local/isofox/main'
include { LINX_REPORT       } from '../modules/local/gpgr/linx_report/main'
include { ORANGE            } from '../modules/local/orange/main'
include { PEACH             } from '../modules/local/peach/main'
include { PROTECT           } from '../modules/local/protect/main'
include { PURPLE            } from '../modules/local/purple/main'
include { SIGS              } from '../modules/local/sigs/main'
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
include { PICARD_COLLECTWGSMETRICS    } from '../modules/nf-core/picard/collectwgsmetrics/main'
include { SAMTOOLS_FLAGSTAT           } from '../modules/nf-core/samtools/flagstat/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Get absolute file paths
samplesheet = Utils.getFileObject(params.input)
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
    if (run.amber || run.cobalt || run.pave || run.lilac) {
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
                def bam = meta.get([Constants.FileType.BAM_WTS, Constants.DataType.TUMOR])
                def meta_isofox = [key: meta.id, id: meta.get(['sample_name', Constants.DataType.TUMOR])]
                return [meta_isofox, bam, "${bam}.bai"]
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
    // MODULE: Run PICARD_COLLECTWGSMETRICS to generate stats required for downstream processes
    //
    if (run.virusinterpreter || run.orange) {
        // Create inputs and create process-specific meta
        // NOTE(SW): CUPPA only requires collectwgsmetrics for the tumor sample in the upstream
        // process Virus Interpreter but ORANGE currently requires collectwgsmetrics for both tumor
        // and normal sample
        // channel: [val(meta_cwm), bam]
        ch_cwm_inputs_all = ch_inputs
            .flatMap { meta ->
                def sample_types = run.orange ? [Constants.DataType.TUMOR, Constants.DataType.NORMAL] : [Constants.DataType.TUMOR]
                return sample_types
                    .collect { sample_type ->
                        def bam = meta.get([Constants.FileType.BAM_WGS, sample_type])
                        def sample_name = meta.get(['sample_name', sample_type])
                        def meta_cwm = [
                            key: meta.id,
                            id: sample_name,
                            // NOTE(SW): must use string representation for caching purposes
                            sample_type_str: sample_type.name(),
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
                def (keys, sample_names, sample_type_strs) = meta_cwm
                    .collect {
                        [it.key, it.id, it.sample_type_str]
                    }
                    .transpose()

                def sample_type_strs_unique = sample_type_strs.unique(false)
                assert sample_type_strs_unique.size() == 1
                def sample_type_str = sample_type_strs_unique[0]

                def meta_cwm_new = [
                    keys: keys,
                    id: sample_names.join('__'),
                    id_simple: keys.join('__'),
                    sample_type_str: sample_type_str,
                ]
                return [meta_cwm_new, *filepaths]
            }

        // Run process
        PICARD_COLLECTWGSMETRICS(
            ch_cwm_inputs,
            PREPARE_REFERENCE.out.genome_fasta,
        )

        // Set outputs, process outputs and restore original meta
        ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS.out.versions)

        // Replicate outputs to reverse unique operation
        // channel: [val(meta_cwm_individual), sample_type_str, metrics]
        ch_cwm_out_individual = PICARD_COLLECTWGSMETRICS.out.metrics
            .flatMap { meta_cwm_shared, metrics ->
                meta_cwm_shared.keys.collect { key ->
                    return [meta_cwm_shared + [key: key], meta_cwm_shared.sample_type_str, metrics]
                }
            }

        // Match outputs to original meta and set output
        // channel (tumor): [val(meta), metrics]
        // channel (normal): [val(meta), metrics]
        ch_cwm_out = WorkflowOncoanalyser.restoreMeta(ch_cwm_out_individual, ch_inputs)
            .branch { meta, sample_type_str, metrics ->
                def sample_type = Utils.getEnumFromString(sample_type_str, Constants.DataType)
                somatic: sample_type == Constants.DataType.TUMOR
                    return [meta, metrics]
                germline: sample_type == Constants.DataType.NORMAL
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
    ch_sage_germline_vcf_out = Channel.empty()
    // channel: [val(meta), sage_coverage]
    ch_sage_germline_coverage_out = Channel.empty()
    // channel: [val(meta), sage_vcf]
    ch_sage_somatic_vcf_out = Channel.empty()
    // channel: [val(meta), bqr_plot]
    ch_sage_somatic_tumor_bqr_out = Channel.empty()
    // channel: [val(meta), bqr_plot]
    ch_sage_somatic_normal_bqr_out = Channel.empty()
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
            hmf_data.sage_coverage_panel_germline,
            hmf_data.sage_highconf_regions,
            hmf_data.sage_pon,
            hmf_data.segment_mappability,
            hmf_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(SAGE.out.versions)
        ch_sage_germline_vcf_out = ch_sage_germline_vcf_out.mix(WorkflowOncoanalyser.restoreMeta(SAGE.out.germline_vcf, ch_inputs))
        ch_sage_germline_coverage_out = ch_sage_germline_coverage_out.mix(WorkflowOncoanalyser.restoreMeta(SAGE.out.germline_coverage, ch_inputs))
        ch_sage_somatic_vcf_out = ch_sage_somatic_vcf_out.mix(WorkflowOncoanalyser.restoreMeta(SAGE.out.somatic_vcf, ch_inputs))
        ch_sage_somatic_tumor_bqr_out = ch_sage_somatic_tumor_bqr_out.mix(WorkflowOncoanalyser.restoreMeta(SAGE.out.somatic_tumor_bqr, ch_inputs))
        ch_sage_somatic_normal_bqr_out = ch_sage_somatic_normal_bqr_out.mix(WorkflowOncoanalyser.restoreMeta(SAGE.out.somatic_normal_bqr, ch_inputs))
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
            ch_pave_germline_inputs_source = ch_sage_germline_vcf_out
            ch_pave_somatic_inputs_source = ch_sage_somatic_vcf_out
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
        ch_purple_inputs_source = WorkflowOncoanalyser.groupByMeta(
            run.amber ? ch_amber_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.AMBER_DIR, Constants.DataType.TUMOR_NORMAL]),
            run.cobalt ? ch_cobalt_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.COBALT_DIR, Constants.DataType.TUMOR_NORMAL]),
            ch_purple_inputs_sv,
            run.pave ? ch_pave_somatic_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.PAVE_VCF, Constants.DataType.TUMOR]),
            run.pave ? ch_pave_germline_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.PAVE_VCF, Constants.DataType.NORMAL]),
        )

        // Create process-specific meta
        // channel: [val(meta_purple), amber_dir, cobalt_dir, sv_hard_vcf, sv_hard_tbi, sv_soft_vcf, sv_soft_tbi, smlv_tumor_vcf, smlv_normal_vcf]
        ch_purple_inputs = ch_purple_inputs_source
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
    // channel: [val(meta), chord_prediction]
    ch_chord_out = Channel.empty()
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
                def meta_chord = [key: meta.id, id: meta.id]
                return [meta_chord, smlv_vcf, sv_vcf]
            }

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
        // Select input sources
        // channel: [val(meta), purple_dir]
        if (run.purple) {
            ch_lilac_inputs_source = ch_purple_out
        } else {
            ch_lilac_inputs_source = WorkflowOncoanalyser.getInput([Constants.FileType.PURPLE_DIR, Constants.DataType.TUMOR_NORMAL])
        }

        // Create channel with all available input BAMs
        // First obtain WTS BAMs
        // channel: [val(meta), wts_bam, wts_bai]
        ch_lilac_bams_wts = Channel.empty()
            .mix(
                ch_inputs_wts_bams = ch_inputs_wts.present.map { meta ->
                    def bam = meta.get([Constants.FileType.BAM_WTS, Constants.DataType.TUMOR])
                    [meta, bam, "${bam}.bai"]
                },
                ch_inputs_wts.absent.map { meta -> [meta, [], []] },
            )

        // Combine WGS and WTS BAMs
        // channel: [val(meta), normal_wgs_bam, normal_wgs_bai, tumor_wgs_bam, tumor_wgs_bai, tumor_wts_bam, tumor_wts_bai]
        ch_lilac_bams = WorkflowOncoanalyser.groupByMeta(
            ch_lilac_bams_wts,
            ch_bams_and_indices,
            flatten_mode: 'nonrecursive',
        )
            .map { data ->
                def meta = data[0]
                def (tbam_wts, tbai_wts, tbam_wgs, nbam_wgs, tbai_wgs, nbai_wgs) = data[1..-1]
                return [meta, nbam_wgs, nbai_wgs, tbam_wgs, tbai_wgs, tbam_wts, tbai_wts]
            }

        // Call subworkflow to run processes
        LILAC(
            ch_lilac_bams,
            run.purple ? ch_purple_out : WorkflowOncoanalyser.getInput([Constants.FileType.PURPLE_DIR, Constants.DataType.TUMOR_NORMAL]),
            PREPARE_REFERENCE.out.genome_fasta,
            PREPARE_REFERENCE.out.genome_fai,
            hmf_data.lilac_resources,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(LILAC.out.versions)
        ch_lilac_out = ch_lilac_out.mix(WorkflowOncoanalyser.restoreMeta(LILAC.out.results, ch_inputs))
    }

    //
    // MODULE: Run VIRUSBreakend and Virus Interpreter to quantify viral content
    //
    ch_virusinterpreter_out = Channel.empty()
    if (run.virusinterpreter) {
        // VIRUSBreakend
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
            ch_cwm_out.somatic,
        )

        // Virus Interpreter
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
        ch_virusinterpreter_out = ch_virusinterpreter_out.mix(WorkflowOncoanalyser.restoreMeta(VIRUSINTERPRETER.out.virusinterpreter, ch_inputs))
    }

    //
    // SUBWORKFLOW: Group structural variants into higher order events with LINX
    //
    // channel: [val(meta), linx_annotation_dir, linx_visuliaser_dir]
    ch_linx_somatic_out = Channel.empty()
    // channel: [val(meta), linx_annotation_dir]
    ch_linx_germline_out = Channel.empty()
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
        ch_linx_germline_out = ch_linx_germline_out.mix(WorkflowOncoanalyser.restoreMeta(LINX.out.germline, ch_inputs))

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
    // MODULE: Run PROTECT to match somatic genomic features with treatment evidence
    //
    // channel: [val(meta), protect]
    ch_protect_out = Channel.empty()
    if (run.protect) {
        // Select input sources
        // channel: [val(meta), linx_somatic_annotation_dir]
        ch_linx_anno = ch_linx_somatic_out.map { meta, anno_dir, vis_dir -> [meta, anno_dir]}

        // channel: [val(meta), chord_prediction, purple_dir, linx_dir, virusinterpreter]
        ch_protect_inputs_source = WorkflowOncoanalyser.groupByMeta(
            run.chord ? ch_chord_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.CHORD_PREDICTION, Constants.DataType.TUMOR]),
            run.purple ? ch_purple_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.PURPLE_DIR, Constants.DataType.TUMOR_NORMAL]),
            run.linx ? ch_linx_anno : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.LINX_ANNO_DIR, Constants.DataType.TUMOR_NORMAL]),
            run.virusinterpreter ? ch_virusinterpreter_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.VIRUSINTERPRETER_TSV, Constants.DataType.TUMOR]),
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
                    tumor_id: meta.get(['sample_name', Constants.DataType.TUMOR]),
                    normal_id: meta.get(['sample_name', Constants.DataType.NORMAL]),
              ]
              return [meta_protect, *other]
            }

        // Run process
        PROTECT(
          ch_protect_inputs,
          PREPARE_REFERENCE.out.genome_version,
          hmf_data.serve_resources,
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
            ch_peach_inputs_source = WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.PURPLE_DIR, Constants.DataType.TUMOR_NORMAL])
        }

        // Create inputs and create process-specific meta
        // channel: [meta_peach, purple_germline_vcf]
        ch_peach_inputs = ch_peach_inputs_source
            .map { meta, purple_dir ->
                def meta_peach = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: meta.get(['sample_name', Constants.DataType.TUMOR]),
                    normal_id: meta.get(['sample_name', Constants.DataType.NORMAL]),
                ]
                def purple_germline_vcf = file(purple_dir).resolve("${meta.get(['sample_name', Constants.DataType.TUMOR])}.purple.germline.vcf.gz")
                file(purple_germline_vcf, checkIfExists: true)
                return [meta_peach, purple_germline_vcf]
            }

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
    // MODULE: Run CUPPA predict tissue of origin
    //
    // channel: [val(meta), cuppa_results]
    ch_cuppa_out = Channel.empty()
    // channel: [val(meta), cuppa_summary_plot]
    ch_cuppa_summary_plot_out = Channel.empty()
    // channel: [val(meta), cuppa_feature_plot]
    ch_cuppa_feature_plot_out = Channel.empty()
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

        // channel: [val(meta), isofox_dir, purple_dir, linx_dir, virusinterpreter]
        ch_cuppa_inputs_source = WorkflowOncoanalyser.groupByMeta(
            ch_cuppa_inputs_isofox,
            run.purple ? ch_purple_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.PURPLE_DIR, Constants.DataType.TUMOR_NORMAL]),
            run.linx ? ch_linx_anno : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.LINX_ANNO_DIR, Constants.DataType.TUMOR_NORMAL]),
            run.virusinterpreter ? ch_virusinterpreter_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.VIRUSINTERPRETER_TSV, Constants.DataType.TUMOR]),
            flatten_mode: 'nonrecursive',
        )

        // Create inputs and create process-specific meta
        // channel: [val(meta_cuppa), isofox_dir, purple_dir, linx_dir, virusinterpreter]
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

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(CUPPA_CLASSIFIER.out.versions)
        ch_cuppa_out = ch_cuppa_out.mix(WorkflowOncoanalyser.restoreMeta(CUPPA_CLASSIFIER.out.csv, ch_inputs))

        // Run process
        CUPPA_VISUALISER(
            CUPPA_CLASSIFIER.out.csv,
        )

        // Set outputs, restoring original meta
        ch_versions = ch_versions.mix(CUPPA_VISUALISER.out.versions)
        ch_cuppa_summary_plot_out = ch_cuppa_summary_plot_out.mix(WorkflowOncoanalyser.restoreMeta(CUPPA_VISUALISER.out.summary_plot, ch_inputs))
        ch_cuppa_feature_plot_out = ch_cuppa_feature_plot_out.mix(WorkflowOncoanalyser.restoreMeta(CUPPA_VISUALISER.out.feature_plot, ch_inputs))
    }

    //
    // MODULE: Run ORANGE to generate static PDF report
    //
    if (run.orange) {
        // SAMtools flagstat
        // Select input source
        // channel (present): [val(meta), sample_type, flagstat]
        // channel (absent): [val(meta)]
        ch_inputs_flagstat = ch_inputs
            .flatMap { meta -> [Constants.DataType.TUMOR, Constants.DataType.NORMAL].collect { [meta, it] } }
            .branch { meta, sample_type ->
                def key = [Constants.FileType.FLAGSTAT, sample_type]
                present: meta.containsKey(key)
                    return [meta, sample_type, meta.get(key)]
                absent: ! meta.containsKey(key)
            }

        // Create inputs and create process-specific meta
        // channel: [val(meta_flagstat), bam, bai]
        ch_flagstat_inputs_all = ch_inputs_flagstat.absent
            .map { meta, sample_type ->
                def bam = meta.get([Constants.FileType.BAM_WGS, sample_type])
                def meta_flagstat = [
                    key: meta.id,
                    id: meta.get(['sample_name', sample_type]),
                    // NOTE(SW): must use string representation for caching purposes
                    sample_type_str: sample_type.name(),
                ]
                return [meta_flagstat, bam, "${bam}.bai"]
            }

        // Collapse duplicate files e.g. repeated normal BAMs for multiple tumor samples
        // channel: [val(meta_flagstat_shared), bam, bai]
        ch_flagstat_inputs = ch_flagstat_inputs_all
            .map { [it[1..-1], it[0]] }
            .groupTuple()
            .map { filepaths, meta_flagstat ->
                def (keys, sample_names, sample_type_strs) = meta_flagstat
                    .collect {
                        [it.key, it.id, it.sample_type_str]
                    }
                    .transpose()

                def sample_type_strs_unique = sample_type_strs.unique(false)
                assert sample_type_strs_unique.size() == 1
                def sample_type_str = sample_type_strs_unique[0]

                def meta_flagstat_new = [
                    keys: keys,
                    id: sample_names.join('__'),
                    id_simple: keys.join('__'),
                    sample_type_str: sample_type_str,
                ]
                return [meta_flagstat_new, *filepaths]
            }

        // Run process
        SAMTOOLS_FLAGSTAT(
            ch_flagstat_inputs,
        )

        // Set version
        ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

        // Replicate outputs to reverse unique operation
        // channel: [val(meta_flagstat_individual), flagstat]
        ch_flagstat_out = SAMTOOLS_FLAGSTAT.out.flagstat
            .flatMap { meta_flagstat_shared, flagstat ->
                def sample_type = Utils.getEnumFromString(meta_flagstat_shared.sample_type_str, Constants.DataType)
                meta_flagstat_shared.keys.collect { key ->
                    return [meta_flagstat_shared + [key: key], sample_type, flagstat]
                }
            }

        // Combine input flagstat channels, restoring original meta where required, split by sample type
        // channel (somatic): [val(meta), flagstat]
        // channel (germline): [val(meta), flagstat]
        ch_orange_inputs_flagstat = Channel.empty()
            .concat(
                ch_inputs_flagstat.present,
                WorkflowOncoanalyser.restoreMeta(ch_flagstat_out, ch_inputs),
            )
            .branch { meta, sample_type, flagstat ->
                somatic: sample_type == Constants.DataType.TUMOR
                    return [meta, flagstat]
                germline: sample_type == Constants.DataType.NORMAL
                    return [meta, flagstat]
            }

        // ORANGE
        // Split LINX channel
        // channel (anno): [val(meta), anno_dir]
        // channel (plot): [val(meta), plot_dir]
        ch_orange_inputs_linx_somatic = ch_linx_somatic_out
            .multiMap { meta, anno_dir, plot_dir ->
                anno: [meta, anno_dir]
                plot: [meta, plot_dir]
            }

        // Select input source
        // NOTE(SW): we could consider not allowing inputs from the samplesheet here since this nothing follows
        ch_orange_inputs_source = WorkflowOncoanalyser.groupByMeta(
            ch_cwm_out.somatic,
            ch_cwm_out.germline,
            ch_orange_inputs_flagstat.somatic,
            ch_orange_inputs_flagstat.germline,
            run.chord ? ch_chord_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.CHORD_PREDICTION, Constants.DataType.TUMOR]),
            run.lilac ? ch_lilac_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.LILAC_DIR, Constants.DataType.NORMAL]),
            run.sage ? ch_sage_somatic_tumor_bqr_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.SAGE_BQR, Constants.DataType.TUMOR]),
            run.sage ? ch_sage_somatic_normal_bqr_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.SAGE_BQR, Constants.DataType.NORMAL]),
            run.sage ? ch_sage_germline_coverage_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.SAGE_COVERAGE, Constants.DataType.NORMAL]),
            run.purple ? ch_purple_out : WorkflowOncoanalyser.getInput([Constants.FileType.PURPLE_DIR, Constants.DataType.TUMOR_NORMAL]),
            run.linx ? ch_orange_inputs_linx_somatic.anno : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.LINX_ANNO_DIR, Constants.DataType.TUMOR_NORMAL]),
            run.linx ? ch_orange_inputs_linx_somatic.plot : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.LINX_PLOT_DIR, Constants.DataType.TUMOR_NORMAL]),
            run.linx ? ch_linx_germline_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.LINX_ANNO_DIR, Constants.DataType.NORMAL]),
            run.protect ? ch_protect_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.PROTECT_TSV, Constants.DataType.TUMOR]),
            run.peach ? ch_peach_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.PEACH_TSV, Constants.DataType.NORMAL]),
            run.cuppa ? ch_cuppa_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.CUPPA_CSV, Constants.DataType.TUMOR]),
            run.cuppa ? ch_cuppa_feature_plot_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.CUPPA_FEATURE_PLOT, Constants.DataType.TUMOR]),
            run.cuppa ? ch_cuppa_summary_plot_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.CUPPA_SUMMARY_PLOT, Constants.DataType.TUMOR]),
            run.virusinterpreter ? ch_virusinterpreter_out : WorkflowOncoanalyser.getInput(ch_inputs, [Constants.FileType.VIRUSINTERPRETER_TSV, Constants.DataType.TUMOR]),
        )

        ch_orange_inputs = ch_orange_inputs_source
            .map {
                def meta = it[0]
                def meta_orange = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: meta.get(['sample_name', Constants.DataType.TUMOR]),
                    normal_id: meta.get(['sample_name', Constants.DataType.NORMAL]),
                ]
                return [meta_orange, *it[1..-1]]
            }

        // Run process
        ORANGE(
            ch_orange_inputs,
            PREPARE_REFERENCE.out.genome_version,
            hmf_data.disease_ontology,
            hmf_data.known_fusion_data,
            hmf_data.driver_gene_panel,
            hmf_data.cohort_mapping,
            hmf_data.cohort_percentiles,
            "5.23 (oncoanalyser 0.0.1dev)",
        )

        // Set outputs
        ch_versions = ch_versions.mix(ORANGE.out.versions)
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
