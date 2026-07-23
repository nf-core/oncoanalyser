//
// PAVE annotates somatic and germline variant VCFs with gene and transcript coding and protein effects
//

include { PAVE_GERMLINE } from '../../../modules/local/pave/germline/main'
include { PAVE_SOMATIC  } from '../../../modules/local/pave/somatic/main'

workflow PAVE_ANNOTATION {
    take:
    // Sample data
    ch_inputs              // channel: [mandatory] [ meta ]
    ch_sage_dir_somatic    // channel: [mandatory] [ meta, sage_dir ]
    ch_sage_dir_germline   // channel: [mandatory] [ meta, sage_dir ]

    // Reference data
    genome_fasta           // channel: [mandatory] /path/to/genome_fasta
    genome_version         // channel: [mandatory] genome version
    genome_fai             // channel: [mandatory] /path/to/genome_fai
    pon_artefacts          // channel: [optional]  /path/to/pon_artefacts
    sage_pon               // channel: [mandatory] /path/to/sage_pon
    sage_blocklist_regions // channel: [mandatory] /path/to/sage_blocklist_regions
    sage_blocklist_sites   // channel: [mandatory] /path/to/sage_blocklist_sites
    clinvar_annotations    // channel: [mandatory] /path/to/clinvar_annotations
    segment_mappability    // channel: [mandatory] /path/to/segment_mappability
    driver_gene_panel      // channel: [mandatory] /path/to/driver_gene_panel
    ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/
    gnomad_resource        // channel: [mandatory] /path/to/gnomad_resource

    // Params
    sequencing_platform    // string:  [mandatory] sequencing platform

    main:
    //
    // MODULE: PAVE germline
    //
    // Select input sources then sort
    // channel: runnable: [ meta, sage_vcf, sage_tbi ]
    // channel: skip: [ meta ]
    ch_sage_germline_inputs_sorted = ch_sage_dir_germline
        .map { meta, sage_dir ->

            def sage_dir_selected = Utils.selectCurrentOrExisting(sage_dir, meta, Constants.INPUT.SAGE_DIR_NORMAL)
            def sage_vcf = sage_dir_selected ? sage_dir_selected.resolve("${Utils.getTumorDnaSampleName(meta)}.sage.germline.vcf.gz") : []
            def sage_tbi = sage_dir_selected ? sage_dir_selected.resolve("${Utils.getTumorDnaSampleName(meta)}.sage.germline.vcf.gz.tbi") : []

            return [meta, sage_vcf, sage_tbi]
        }
        .branch { meta, sage_vcf, sage_tbi ->

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.PAVE_DIR_NORMAL)

            runnable: Utils.hasTumorDna(meta) && Utils.hasNormalDna(meta) && sage_vcf && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_pave, sage_vcf, sage_tbi ]
    ch_pave_germline_inputs = ch_sage_germline_inputs_sorted.runnable
        .map { meta, sage_vcf, sage_tbi ->

            def meta_pave = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
            ]

            return [meta_pave, sage_vcf, sage_tbi]
        }

    // Run process
    PAVE_GERMLINE(
        ch_pave_germline_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        sage_blocklist_regions,
        sage_blocklist_sites,
        clinvar_annotations,
        segment_mappability,
        driver_gene_panel,
        ensembl_data_resources,
        sequencing_platform,
    )

    //
    // MODULE: PAVE somatic
    //
    // Select input sources then sort
    // channel: runnable: [ meta, sage_vcf, sage_tbi ]
    // channel: skip: [ meta ]
    ch_sage_somatic_inputs_sorted = ch_sage_dir_somatic
        .map { meta, sage_dir ->

            def sage_dir_selected = Utils.selectCurrentOrExisting(sage_dir, meta, Constants.INPUT.SAGE_DIR_TUMOR)
            def sage_vcf = sage_dir_selected ? sage_dir_selected.resolve("${Utils.getTumorDnaSampleName(meta)}.sage.somatic.vcf.gz") : []
            def sage_tbi = sage_dir_selected ? sage_dir_selected.resolve("${Utils.getTumorDnaSampleName(meta)}.sage.somatic.vcf.gz.tbi") : []

            return [meta, sage_vcf, sage_tbi]
        }
        .branch { meta, sage_vcf, sage_tbi ->

            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.PAVE_DIR_TUMOR)

            runnable: Utils.hasTumorDna(meta) && sage_vcf && ! has_existing
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_pave, sage_vcf, sage_tbi ]
    ch_pave_somatic_inputs = ch_sage_somatic_inputs_sorted.runnable
        .map { meta, sage_vcf, sage_tbi ->

            def meta_pave = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorDnaSampleName(meta),
            ]

            return [meta_pave, sage_vcf, sage_tbi]
        }

    // Run process
    PAVE_SOMATIC(
        ch_pave_somatic_inputs,
        genome_fasta,
        genome_version,
        genome_fai,
        pon_artefacts,
        sage_pon,
        clinvar_annotations,
        segment_mappability,
        driver_gene_panel,
        ensembl_data_resources,
        gnomad_resource,
        sequencing_platform,
    )

    //
    // STEP: Outputs
    //
    // Set outputs, restoring original meta
    // channel: [ meta, pave_dir ]
    ch_outputs_somatic = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('pave_somatic_dir'), ch_inputs),
            ch_sage_somatic_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    // channel: [ meta, pave_dir ]
    ch_outputs_germline = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('pave_germline_dir'), ch_inputs),
            ch_sage_germline_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    germline_dir = ch_outputs_germline // channel: [ meta, pave_dir ]
    somatic_dir  = ch_outputs_somatic  // channel: [ meta, pave_dir ]
}
