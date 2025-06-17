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
    params.isofox_counts,
    params.isofox_gc_ratios,
    params.isofox_gene_ids,
    params.isofox_tpm_norm,
]

if (run_config.stages.lilac) {
    if (params.genome_version.toString() == '38' && params.genome_type == 'alt' && params.containsKey('ref_data_hla_slice_bed')) {
        checkPathParamList.add(params.ref_data_hla_slice_bed)
    }
}

// TODO(SW): consider whether we should check for null entries here for errors to be more informative
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Used in Isofox subworkflow only
isofox_read_length = params.isofox_read_length !== null ? params.isofox_read_length : Constants.DEFAULT_ISOFOX_READ_LENGTH_TARGETED

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { AMBER_PROFILING       } from '../subworkflows/local/amber_profiling'
include { BAMTOOLS_METRICS      } from '../subworkflows/local/bamtools_metrics'
include { CIDER_CALLING         } from '../subworkflows/local/cider_calling'
include { COBALT_PROFILING      } from '../subworkflows/local/cobalt_profiling'
include { ESVEE_CALLING         } from '../subworkflows/local/esvee_calling'
include { ISOFOX_QUANTIFICATION } from '../subworkflows/local/isofox_quantification'
include { LILAC_CALLING         } from '../subworkflows/local/lilac_calling'
include { LINX_ANNOTATION       } from '../subworkflows/local/linx_annotation'
include { LINX_PLOTTING         } from '../subworkflows/local/linx_plotting'
include { ORANGE_REPORTING      } from '../subworkflows/local/orange_reporting'
include { PAVE_ANNOTATION       } from '../subworkflows/local/pave_annotation'
include { PEACH_CALLING         } from '../subworkflows/local/peach_calling'
include { PREPARE_REFERENCE     } from '../subworkflows/local/prepare_reference'
include { PURPLE_CALLING        } from '../subworkflows/local/purple_calling'
include { READ_ALIGNMENT_DNA    } from '../subworkflows/local/read_alignment_dna'
include { READ_ALIGNMENT_RNA    } from '../subworkflows/local/read_alignment_rna'
include { REDUX_PROCESSING      } from '../subworkflows/local/redux_processing'
include { SAGE_APPEND           } from '../subworkflows/local/sage_append'
include { SAGE_CALLING          } from '../subworkflows/local/sage_calling'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Get absolute file paths
samplesheet = Utils.getFileObject(params.input)

workflow TARGETED {
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
    panel_data = PREPARE_REFERENCE.out.panel_data

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

        READ_ALIGNMENT_RNA(
            ch_inputs,
            ref_data.genome_star_index,
        )

        ch_versions = ch_versions.mix(
            READ_ALIGNMENT_DNA.out.versions,
            READ_ALIGNMENT_RNA.out.versions,
        )

        ch_align_dna_tumor_out = ch_align_dna_tumor_out.mix(READ_ALIGNMENT_DNA.out.dna_tumor)
        ch_align_dna_normal_out = ch_align_dna_normal_out.mix(READ_ALIGNMENT_DNA.out.dna_normal)
        ch_align_dna_donor_out = ch_align_dna_donor_out.mix(READ_ALIGNMENT_DNA.out.dna_donor)
        ch_align_rna_tumor_out = ch_align_rna_tumor_out.mix(READ_ALIGNMENT_RNA.out.rna_tumor)

    } else {

        ch_align_dna_tumor_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_align_dna_normal_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_align_dna_donor_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_align_rna_tumor_out = ch_inputs.map { meta -> [meta, [], []] }

    }

    //
    // SUBWORKFLOW: Run REDUX for DNA BAMs
    //
    // channel: [ meta, bam, bai ]
    ch_redux_dna_tumor_out = Channel.empty()
    ch_redux_dna_normal_out = Channel.empty()
    ch_redux_dna_donor_out = Channel.empty()

    // channel: [ meta, dup_freq_tsv, jitter_tsv, ms_tsv ]
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
    // MODULE: Run Isofox to analyse RNA data
    //
    // channel: [ meta, isofox_dir ]
    ch_isofox_out = Channel.empty()
    if (run_config.stages.isofox) {

        isofox_counts = params.isofox_counts ? file(params.isofox_counts) : panel_data.isofox_counts
        isofox_gc_ratios = params.isofox_gc_ratios ? file(params.isofox_gc_ratios) : panel_data.isofox_gc_ratios

        isofox_gene_ids = params.isofox_gene_ids ? file(params.isofox_gene_ids) : panel_data.isofox_gene_ids
        isofox_tpm_norm = params.isofox_tpm_norm ? file(params.isofox_tpm_norm) : panel_data.isofox_tpm_norm

        ISOFOX_QUANTIFICATION(
            ch_inputs,
            ch_align_rna_tumor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            hmf_data.ensembl_data_resources,
            hmf_data.known_fusion_data,
            isofox_counts,
            isofox_gc_ratios,
            isofox_gene_ids,
            isofox_tpm_norm,
            params.isofox_functions,
            isofox_read_length,
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
            ch_redux_dna_tumor_out,
            ch_redux_dna_normal_out,
            ch_redux_dna_donor_out,
            ref_data.genome_version,
            hmf_data.heterozygous_sites,
            panel_data.target_region_bed,
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
            panel_data.target_region_normalisation,
        )

        ch_versions = ch_versions.mix(COBALT_PROFILING.out.versions)

        ch_cobalt_out = ch_cobalt_out.mix(COBALT_PROFILING.out.cobalt_dir)

    } else {

        ch_cobalt_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Call structural variants with ESVEE
    //
    // channel: [ meta, esvee_vcf ]
    ch_esvee_germline_out = Channel.empty()
    ch_esvee_somatic_out = Channel.empty()
    if (run_config.stages.esvee) {

        ESVEE_CALLING(
            ch_inputs,
            ch_redux_dna_tumor_out,
            ch_redux_dna_normal_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            ref_data.genome_img,
            hmf_data.sv_prep_blocklist,
            hmf_data.known_fusions,
            hmf_data.gridss_pon_breakends,
            hmf_data.gridss_pon_breakpoints,
            hmf_data.repeatmasker_annotations,
            hmf_data.decoy_sequences_image,
            hmf_data.unmap_regions,
        )

        ch_versions = ch_versions.mix(ESVEE_CALLING.out.versions)

        ch_esvee_germline_out = ch_esvee_germline_out.mix(ESVEE_CALLING.out.germline_vcf)
        ch_esvee_somatic_out = ch_esvee_somatic_out.mix(ESVEE_CALLING.out.somatic_vcf)

    } else {

        ch_esvee_germline_out = ch_inputs.map { meta -> [meta, []] }
        ch_esvee_somatic_out = ch_inputs.map { meta -> [meta, []] }

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
            ch_redux_dna_tumor_out,
            ch_redux_dna_normal_out,
            ch_redux_dna_donor_out,
            ch_redux_dna_tumor_tsv_out,
            ch_redux_dna_normal_tsv_out,
            ch_redux_dna_donor_tsv_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            hmf_data.sage_known_hotspots_somatic,
            hmf_data.sage_known_hotspots_germline,
            panel_data.sage_actionable_panel,
            panel_data.sage_coverage_panel,
            hmf_data.sage_highconf_regions,
            hmf_data.segment_mappability,
            panel_data.driver_gene_panel,
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
            panel_data.pon_artefacts,
            hmf_data.sage_blocklist_regions,
            hmf_data.sage_blocklist_sites,
            hmf_data.clinvar_annotations,
            hmf_data.segment_mappability,
            panel_data.driver_gene_panel,
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
            ch_esvee_somatic_out,
            ch_esvee_germline_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            hmf_data.gc_profile,
            hmf_data.sage_known_hotspots_somatic,
            hmf_data.sage_known_hotspots_germline,
            panel_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
            hmf_data.purple_germline_del,
            panel_data.target_region_bed,
            panel_data.target_region_ratios,
            panel_data.target_region_msi_indels,
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
            ch_align_rna_tumor_out,
            ch_purple_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            true,  // run_germline
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
            panel_data.driver_gene_panel,
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
    // channel: [ meta, linx_visualiser_dir ]
    ch_linx_somatic_visualiser_dir_out = Channel.empty()
    if (run_config.stages.linx) {

        LINX_PLOTTING(
            ch_inputs,
            ch_linx_somatic_out,
            ref_data.genome_version,
            hmf_data.ensembl_data_resources,
        )

        ch_versions = ch_versions.mix(LINX_PLOTTING.out.versions)

        ch_linx_somatic_visualiser_dir_out = ch_linx_somatic_visualiser_dir_out.mix(LINX_PLOTTING.out.visualiser_dir)

    } else {

        ch_linx_somatic_visualiser_dir_out = ch_inputs.map { meta -> [meta, []] }

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
            ch_redux_dna_tumor_out,
            ch_redux_dna_normal_out,
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
    // SUBWORKFLOW: Run CIDER to identify and annotate CDR3 sequences of IG and TCR loci
    //
    if (run_config.stages.cider) {

        CIDER_CALLING(
            ch_inputs,
            ch_redux_dna_tumor_out,
            ch_align_rna_tumor_out,
            ref_data.genome_version,
            hmf_data.cider_blastdb,
        )

        ch_versions = ch_versions.mix(CIDER_CALLING.out.versions)

    }

    //
    // SUBWORKFLOW: Run LILAC for HLA typing and somatic CNV and SNV calling
    //
    // channel: [ meta, lilac_dir ]
    ch_lilac_out = Channel.empty()
    if (run_config.stages.lilac) {

        // Set HLA slice BED if provided in params
        ref_data_hla_slice_bed = params.containsKey('ref_data_hla_slice_bed') ? params.ref_data_hla_slice_bed : []

        LILAC_CALLING(
            ch_inputs,
            ch_redux_dna_tumor_out,
            ch_redux_dna_normal_out,
            ch_align_rna_tumor_out,
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
    // SUBWORKFLOW: Run PEACH to call germline haplotypes and report pharmacogenomics
    //
    // channel: [ meta, peach_dir ]
    ch_peach_out = Channel.empty()
    if (run_config.stages.peach) {

        PEACH_CALLING(
            ch_inputs,
            ch_purple_out,
            hmf_data.peach_haplotypes,
            hmf_data.peach_haplotype_functions,
            hmf_data.peach_drug_info,
        )

        ch_versions = ch_versions.mix(PEACH_CALLING.out.versions)

        ch_peach_out = ch_peach_out.mix(PEACH_CALLING.out.peach_dir)

    } else {

        ch_peach_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run ORANGE to generate static PDF report
    //
    if (run_config.stages.orange) {

        // Create placeholder channels for empty remaining channels
        ch_chord_out = ch_inputs.map { meta -> [meta, []] }
        ch_cuppa_out = ch_inputs.map { meta -> [meta, []] }
        ch_sigs_out = ch_inputs.map { meta -> [meta, []] }
        ch_virusinterpreter_out = ch_inputs.map { meta -> [meta, []] }

        ORANGE_REPORTING(
            ch_inputs,
            ch_bamtools_somatic_out,
            ch_bamtools_germline_out,
            ch_sage_somatic_dir_out,
            ch_sage_germline_dir_out,
            ch_sage_somatic_append_out,
            ch_sage_germline_append_out,
            ch_purple_out,
            ch_linx_somatic_out,
            ch_linx_somatic_visualiser_dir_out,
            ch_linx_germline_out,
            ch_virusinterpreter_out,
            ch_chord_out,
            ch_sigs_out,
            ch_lilac_out,
            ch_cuppa_out,
            ch_peach_out,
            ch_isofox_out,
            ref_data.genome_version,
            hmf_data.disease_ontology,
            hmf_data.cohort_mapping,
            hmf_data.cohort_percentiles,
            hmf_data.known_fusion_data,
            panel_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
            hmf_data.sigs_etiology,
            hmf_data.alt_sj_distribution,
            hmf_data.gene_exp_distribution,
        )

        ch_versions = ch_versions.mix(ORANGE_REPORTING.out.versions)
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
