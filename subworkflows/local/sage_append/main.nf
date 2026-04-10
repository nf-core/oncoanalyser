//
// SAGE append adds additional sample data to an existing SAGE VCF
//

include { SAGE_APPEND as SAGE_APPEND_SOMATIC  } from '../../../modules/local/sage/append/main'
include { SAGE_APPEND as SAGE_APPEND_GERMLINE } from '../../../modules/local/sage/append/main'

workflow SAGE_APPEND {
    take:
    // Sample data
    ch_inputs         // channel: [mandatory] [ meta ]
    ch_purple_dir     // channel: [mandatory] [ meta, purple_dir ]
    ch_tumor_dna_bam  // channel: [mandatory] [ meta, bam, bai ]
    ch_tumor_dna_tsv  // channel: [mandatory] [ meta, redux_tsv, ... ]
    ch_tumor_rna_bam  // channel: [mandatory] [ meta, bam, bai ]

    // Reference data
    genome_fasta     // channel: [mandatory] /path/to/genome_fasta
    genome_version   // channel: [mandatory] genome version
    genome_fai       // channel: [mandatory] /path/to/genome_fai
    genome_dict      // channel: [mandatory] /path/to/genome_dict

    // Params
    sequencing_type  // string:  [mandatory] sequencing type
    enable_germline  // boolean: [mandatory] Enable germline
    targeted_mode    // boolean: [mandatory] Set targeted mode

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    def run_mode = pipeline.PipelineMode.fromString(params.mode)
    def purity_estimate_mode = run_mode === pipeline.PipelineMode.PURITY_ESTIMATE

    // Select input sources and sort
    // channel: runnable: { meta, tumor_dna_bam, tumor_dna_bai, [tumor_dna_redux_tsv, ...], tumor_rna_bam, tumor_rna_bai, purple_dir }
    // channel: skip: [ meta ]
    ch_inputs_sorted = channels.WorkflowChannels.groupByMeta(
        flatten_mode: 'singletons_only',
        ch_tumor_dna_bam,
        ch_tumor_dna_tsv,
        ch_tumor_rna_bam,
        ch_purple_dir,
    )
        .map { meta, tumor_dna_bam_bai, tumor_dna_tsvs, tumor_rna_bam_bai, purple_dir ->

            def (tumor_dna_bam, tumor_dna_bai) = Inputs.resolveReduxBamBai(tumor_dna_bam_bai, meta, samplesheet.SampleType.TUMOR)
            def tumor_dna_redux_tsvs = Inputs.resolveReduxTsvFiles(tumor_dna_tsvs, meta, samplesheet.SampleType.TUMOR)

            def (tumor_rna_bam, tumor_rna_bai) = tumor_rna_bam_bai
            tumor_rna_bam = Inputs.preferUserProvidedInput(tumor_rna_bam, meta, Inputs.KEY.BAM_RNA_TUMOR)
            tumor_rna_bai = Inputs.preferPipelineOutput(tumor_rna_bai, meta, Inputs.KEY.BAI_RNA_TUMOR)

            purple_dir = Inputs.preferUserProvidedInput(purple_dir, meta, Inputs.KEY.PURPLE_DIR)

            def inputs = [
                meta: meta,
                tumor_dna_bam: tumor_dna_bam,
                tumor_dna_bai: tumor_dna_bai,
                tumor_dna_redux_tsvs: tumor_dna_redux_tsvs,
                tumor_rna_bam: tumor_rna_bam,
                tumor_rna_bai: tumor_rna_bai,
                purple_dir: purple_dir
            ]
        }
        .branch { inputs ->
            def has_bam = inputs.tumor_dna_bam || inputs.tumor_rna_bam
            runnable: has_bam && inputs.purple_dir
                return inputs
            skip: true
                return inputs.meta
        }

    //
    // MODULE: SAGE append germline
    //
    // Select inputs that are eligible to run
    // channel: runnable: { meta, tumor_dna_bam, tumor_dna_bai, [tumor_dna_redux_tsv, ...], tumor_rna_bam, tumor_rna_bai, purple_dir }
    // channel: skip: [ meta ]
    ch_inputs_germline_sorted = ch_inputs_sorted.runnable
        .branch { inputs ->

            def meta = inputs.meta

            def has_tumor_rna = Inputs.hasTumorRna(meta)
            def has_normal_dna = Inputs.hasNormalDna(meta)
            def has_smlv_germline = Inputs.resolvePurpleGermlineVcf(inputs.purple_dir, meta)

            def should_append_rna_variants = has_tumor_rna && has_normal_dna && has_smlv_germline

            def has_existing = Inputs.hasExistingInput(meta, Inputs.KEY.SAGE_APPEND_DIR_NORMAL)

            runnable: should_append_rna_variants && !has_existing && enable_germline
                return inputs
            skip: true
                return inputs.meta
        }

    // Create process input channel
    // channel: [ meta_append, purple_smlv_vcf, [bam, ...], [bai, ...], [tumor_dna_redux_tsv, ...] ]
    ch_sage_append_germline_inputs = ch_inputs_germline_sorted.runnable
        .map { inputs ->

            def meta = inputs.meta

            // NOTE(SW): explicit in expectation to always obtain the primary tumor DNA sample ID here
            def tumor_dna_id = Inputs.getTumorDnaSampleName(meta, 'primary')
            def output_file_id = Inputs.getNormalDnaSampleName(meta)

            def meta_append = [
                key: meta.group_id,
                id: meta.group_id,
                output_file_id: output_file_id,
                reference_ids: [Inputs.getTumorRnaSampleName(meta)],
            ]

            def bams = [inputs.tumor_rna_bam]
            def bais = [inputs.tumor_rna_bai]

            def purple_smlv_vcf = Inputs.resolvePurpleGermlineVcf(inputs.purple_dir, meta)

            return [meta_append, purple_smlv_vcf, bams, bais, []]
        }

    // Run process
    SAGE_APPEND_GERMLINE(
        ch_sage_append_germline_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        sequencing_type,
        targeted_mode,
    )

    ch_versions = ch_versions.mix(SAGE_APPEND_GERMLINE.out.versions)

    //
    // MODULE: SAGE append somatic
    //
    // Select inputs that are eligible to run
    // channel: runnable: { meta, tumor_dna_bam, tumor_dna_bai, [tumor_dna_redux_tsv, ...], tumor_rna_bam, tumor_rna_bai, purple_dir }
    // channel: skip: [ meta ]
    ch_inputs_somatic_sorted = ch_inputs_sorted.runnable
        .branch { inputs ->

            def meta = inputs.meta

            def has_tumor_rna = Inputs.hasTumorRna(meta)
            def has_tumor_dna = Inputs.hasTumorDna(meta)
            def has_smlv_somatic = Inputs.resolvePurpleSomaticVcf(inputs.purple_dir, meta)

            def should_append_rna_variants = !purity_estimate_mode && has_tumor_rna && has_tumor_dna && has_smlv_somatic
            def should_append_longitudinal_variants = purity_estimate_mode && has_tumor_dna && has_smlv_somatic

            def has_existing = Inputs.hasExistingInput(meta, Inputs.KEY.SAGE_APPEND_DIR_TUMOR)

            runnable: (should_append_rna_variants || should_append_longitudinal_variants) && !has_existing
                return inputs
            skip: true
                return inputs.meta
        }

    // Create process input channel
    // channel: [ meta_append, purple_smlv_vcf, [bam, ...], [bai, ...], [tumor_dna_redux_tsv, ...] ]
    ch_sage_append_somatic_inputs = ch_inputs_somatic_sorted.runnable
        .map { inputs ->

            def meta = inputs.meta

            def output_file_id = purity_estimate_mode
                ? Inputs.getTumorDnaSampleName(meta, 'longitudinal')
                : Inputs.getTumorDnaSampleName(meta, 'primary')

            def meta_append = [
                key: meta.group_id,
                id: meta.group_id,
                output_file_id: output_file_id,
                reference_ids: [],
            ]

            def bams = []
            def bais = []
            def redux_tsvs = []

            if (!purity_estimate_mode && inputs.tumor_rna_bam) {
                meta_append.reference_ids.add(Inputs.getTumorRnaSampleName(meta))
                bams.add(inputs.tumor_rna_bam)
                bais.add(inputs.tumor_rna_bai)
            }

            if (purity_estimate_mode && inputs.tumor_dna_bam) {
                meta_append.reference_ids.add(Inputs.getTumorDnaSampleName(meta, 'longitudinal'))
                bams.add(inputs.tumor_dna_bam)
                bais.add(inputs.tumor_dna_bai)
                redux_tsvs = inputs.tumor_dna_redux_tsvs
            }

            def purple_smlv_vcf = Inputs.resolvePurpleSomaticVcf(inputs.purple_dir, meta, 'primary')

            return [meta_append, purple_smlv_vcf, bams, bais, redux_tsvs]
        }

    // Run process
    SAGE_APPEND_SOMATIC(
        ch_sage_append_somatic_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        genome_dict,
        sequencing_type,
        targeted_mode,
    )

    ch_versions = ch_versions.mix(SAGE_APPEND_SOMATIC.out.versions)

    // Set outputs, restoring original meta
    // channel: [ meta, sage_append_dir ]
    ch_somatic_dir = Channel.empty()
        .mix(
            channels.WorkflowChannels.restoreMeta(SAGE_APPEND_SOMATIC.out.sage_append_dir, ch_inputs),
            channels.PlaceholderChannels.toolDir(ch_inputs_somatic_sorted.skip),
            channels.PlaceholderChannels.toolDir(ch_inputs_sorted.skip),
        )

    ch_germline_dir = Channel.empty()
        .mix(
            channels.WorkflowChannels.restoreMeta(SAGE_APPEND_GERMLINE.out.sage_append_dir, ch_inputs),
            channels.PlaceholderChannels.toolDir(ch_inputs_germline_sorted.skip),
            channels.PlaceholderChannels.toolDir(ch_inputs_sorted.skip),
        )

    emit:
    somatic_dir  = ch_somatic_dir  // channel: [ meta, sage_append_dir ]
    germline_dir = ch_germline_dir // channel: [ meta, sage_append_dir ]

    versions     = ch_versions     // channel: [ versions.yml ]
}
