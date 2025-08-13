import Constants
import Processes
import Utils

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AMBER_PROFILING            } from '../subworkflows/local/amber_profiling'
include { COBALT_NORMALISATION       } from '../subworkflows/local/cobalt_normalisation'
include { COBALT_PROFILING           } from '../subworkflows/local/cobalt_profiling'
include { ISOFOX_NORMALISATION       } from '../subworkflows/local/isofox_normalisation'
include { ISOFOX_QUANTIFICATION      } from '../subworkflows/local/isofox_quantification'
include { PAVE_PON_CREATION          } from '../subworkflows/local/pave_pon_creation'
include { PREPARE_REFERENCE          } from '../subworkflows/local/prepare_reference'
include { READ_ALIGNMENT_DNA         } from '../subworkflows/local/read_alignment_dna'
include { READ_ALIGNMENT_RNA         } from '../subworkflows/local/read_alignment_rna'
include { REDUX_PROCESSING           } from '../subworkflows/local/redux_processing'
include { SAGE_CALLING               } from '../subworkflows/local/sage_calling'

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PANEL_RESOURCE_CREATION {
    take:
    inputs
    run_config

    main:
    // Check input path parameters to see if they exist
    def checkPathParamList = [
        params.isofox_counts,
        params.isofox_gc_ratios,
        params.isofox_gene_ids,
        params.isofox_tpm_norm,
        params.driver_gene_panel,
        params.target_regions_bed,
    ]

    if (run_config.stages.lilac) {
        if (params.genome_version.toString() == '38' && params.genome_type == 'alt' && params.containsKey('ref_data_hla_slice_bed')) {
            checkPathParamList.add(params.ref_data_hla_slice_bed)
        }
    }

    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

    // Create channel for versions
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Create input channel from parsed CSV
    // channel: [ meta ]
    ch_inputs = Channel.fromList(inputs)

    // Set up reference data, assign more human readable variables
    prep_config = WorkflowMain.getPrepConfigForRun(run_config)
    PREPARE_REFERENCE(
        prep_config,
    )
    ref_data = PREPARE_REFERENCE.out
    hmf_data = PREPARE_REFERENCE.out.hmf_data

    ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)

    //
    // SUBWORKFLOW: Run read alignment to generate BAMs
    //
    READ_ALIGNMENT_DNA(
        ch_inputs,
        ref_data.genome_fasta,
        ref_data.genome_bwamem2_index,
        params.max_fastq_records,
        params.fastp_umi_enabled,
        params.fastp_umi_location,
        params.fastp_umi_length,
        params.fastp_umi_skip,
    )

    READ_ALIGNMENT_RNA(
        ch_inputs,
        ref_data.genome_star_index,
    )

    // channel: [ meta, [bam, ...], [bai, ...] ]
    ch_versions = ch_versions.mix(
        READ_ALIGNMENT_DNA.out.versions,
        READ_ALIGNMENT_RNA.out.versions,
    )

    // channel: [ meta, [bam, ...], [bai, ...] ]
    ch_align_dna_tumor_out = READ_ALIGNMENT_DNA.out.dna_tumor
    ch_align_dna_normal_out = READ_ALIGNMENT_DNA.out.dna_normal
    ch_align_rna_tumor_out = READ_ALIGNMENT_RNA.out.rna_tumor

    //
    // SUBWORKFLOW: Run REDUX for DNA BAMs
    //
    REDUX_PROCESSING(
        ch_inputs,
        ch_align_dna_tumor_out,
        ch_align_dna_normal_out,
        ch_inputs.map { meta -> [meta, [], []] },  // ch_dna_donor
        ref_data.genome_fasta,
        ref_data.genome_version,
        ref_data.genome_fai,
        ref_data.genome_dict,
        hmf_data.unmap_regions,
        hmf_data.msi_jitter_sites,
        params.redux_umi_enabled,
        params.redux_umi_duplex_delim,
    )

    ch_versions = ch_versions.mix(REDUX_PROCESSING.out.versions)

    // channel: [ meta, bam, bai ]
    ch_redux_dna_tumor_out = REDUX_PROCESSING.out.dna_tumor
    ch_redux_dna_normal_out = REDUX_PROCESSING.out.dna_normal

    // channel: [ meta, dup_freq_tsv, jitter_tsv, ms_tsv, repeat_tsv ]
    ch_redux_dna_tumor_tsv_out = REDUX_PROCESSING.out.dna_tumor_tsv
    ch_redux_dna_normal_tsv_out = REDUX_PROCESSING.out.dna_normal_tsv

    //
    // MODULE: Run Isofox to analyse RNA data
    //
    isofox_counts = params.isofox_counts ? file(params.isofox_counts) : hmf_data.isofox_counts
    isofox_gc_ratios = params.isofox_gc_ratios ? file(params.isofox_gc_ratios) : hmf_data.isofox_gc_ratios
    isofox_read_length = params.isofox_read_length !== null ? params.isofox_read_length : Constants.DEFAULT_ISOFOX_READ_LENGTH_TARGETED

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
        [],  // isofox_gene_ids
        [],  // isofox_tpm_norm
        'TRANSCRIPT_COUNTS',
        isofox_read_length,
    )

    ch_versions = ch_versions.mix(ISOFOX_QUANTIFICATION.out.versions)

    // channel: [ meta, isofox_dir ]
    ch_isofox_out = ISOFOX_QUANTIFICATION.out.isofox_dir

    //
    // SUBWORKFLOW: Run AMBER to obtain b-allele frequencies
    //
    AMBER_PROFILING(
        ch_inputs,
        ch_redux_dna_tumor_out,
        ch_redux_dna_normal_out,
        ch_inputs.map { meta -> [meta, [], []] },  // ch_donor_bam
        ref_data.genome_version,
        hmf_data.heterozygous_sites,
        params.target_regions_bed,
        2,   // tumor_min_depth
    )

    ch_versions = ch_versions.mix(AMBER_PROFILING.out.versions)

    // channel: [ meta, amber_dir ]
    ch_amber_out = AMBER_PROFILING.out.amber_dir

    //
    // SUBWORKFLOW: Run COBALT to obtain read ratios
    //
    COBALT_PROFILING(
        ch_inputs,
        ch_redux_dna_tumor_out,
        ch_redux_dna_normal_out,
        hmf_data.gc_profile,
        hmf_data.diploid_bed,
        [],  // panel_target_region_normalisation
        true,  // targeted_mode
    )

    ch_versions = ch_versions.mix(COBALT_PROFILING.out.versions)

    // channel: [ meta, cobalt_dir ]
    ch_cobalt_out = COBALT_PROFILING.out.cobalt_dir


    // SUBWORKFLOW: call SNV, MNV, and small INDELS with SAGE
    //
    SAGE_CALLING(
        ch_inputs,
        ch_redux_dna_tumor_out,
        ch_redux_dna_normal_out,
        ch_inputs.map { meta -> [meta, [], []] },  // ch_donor_bam
        ch_redux_dna_tumor_tsv_out,
        ch_redux_dna_normal_tsv_out,
        ch_inputs.map { meta -> [meta, [], [], []] },  // ch_donor_tsv
        ref_data.genome_fasta,
        ref_data.genome_version,
        ref_data.genome_fai,
        ref_data.genome_dict,
        hmf_data.sage_pon,
        hmf_data.sage_known_hotspots_somatic,
        hmf_data.sage_known_hotspots_germline,
        hmf_data.sage_highconf_regions,
        hmf_data.segment_mappability,
        params.driver_gene_panel,
        hmf_data.ensembl_data_resources,
        hmf_data.gnomad_resource,
        true,  // enable_germline
        true,  // targeted_mode
    )

    ch_versions = ch_versions.mix(SAGE_CALLING.out.versions)

    // channel: [ meta, sage_vcf, sage_tbi ]
    ch_sage_somatic_vcf_out = SAGE_CALLING.out.somatic_vcf

    //
    // SUBWORKFLOW: Run COBALT normalisation
    //
    target_regions_bed = params.target_regions_bed ? file(params.target_regions_bed) : []

    COBALT_NORMALISATION(
        ch_amber_out,
        ch_cobalt_out,
        ref_data.genome_version,
        hmf_data.gc_profile,
        target_regions_bed,
    )

    ch_versions = ch_versions.mix(COBALT_NORMALISATION.out.versions)

    //
    // SUBWORKFLOW: Run PAVE panel of normals creation
    //
    PAVE_PON_CREATION(
        ch_sage_somatic_vcf_out,
        ref_data.genome_version,
    )

    ch_versions = ch_versions.mix(PAVE_PON_CREATION.out.versions)

    //
    // SUBWORKFLOW: Run Isofox TPM normalisation
    //
    isofox_gene_ids = params.isofox_gene_ids ? file(params.isofox_gene_ids) : []

    ISOFOX_NORMALISATION(
        ch_isofox_out,
        ref_data.genome_version,
        isofox_gene_ids,
        hmf_data.gene_exp_distribution,
    )

    ch_versions = ch_versions.mix(ISOFOX_NORMALISATION.out.versions)

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
