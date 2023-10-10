//
// PAVE annotates somatic and germline variant VCFs with gene and transcript coding and protein effects
//

import Constants
import Utils

include { PAVE_GERMLINE as GERMLINE } from '../../modules/local/pave/germline/main'
include { PAVE_SOMATIC as SOMATIC   } from '../../modules/local/pave/somatic/main'

workflow PAVE_ANNOTATION {
    take:
        // Sample data
        ch_inputs              // channel: [mandatory] [ meta ]
        ch_sage_germline_vcf   // channel: [mandatory] [ meta, sage_germline_vcf, sage_somatic_tbi ]
        ch_sage_somatic_vcf    // channel: [mandatory] [ meta, sage_somatic_vcf, sage_somatic_tbi ]

        // Reference data
        genome_fasta           // channel: [mandatory] /path/to/genome_fasta
        genome_version         // channel: [mandatory] genome version
        genome_fai             // channel: [mandatory] /path/to/genome_fai
        sage_pon               // channel: [mandatory] /path/to/sage_pon
        pon_artefacts          // channel: [optional]  /path/to/pon_artefacts
        sage_blocklist_regions // channel: [mandatory] /path/to/sage_blocklist_regions
        sage_blocklist_sites   // channel: [mandatory] /path/to/sage_blocklist_sites
        clinvar_annotations    // channel: [mandatory] /path/to/clinvar_annotations
        segment_mappability    // channel: [mandatory] /path/to/segment_mappability
        driver_gene_panel      // channel: [mandatory] /path/to/driver_gene_panel
        ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/
        gnomad_resource        // channel: [mandatory] /path/to/gnomad_resource

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        //
        // MODULE: PAVE germline
        //
        // Select input sources and sort
        // channel: runnable: [ meta, sage_vcf ]
        // channel: skip: [ meta ]
        ch_sage_germline_inputs_sorted = ch_sage_germline_vcf
            .map { meta, sage_vcf, sage_tbi ->
                return [
                    meta,
                    Utils.selectCurrentOrExisting(sage_vcf, meta, Constants.INPUT.SAGE_VCF_NORMAL),
                ]
            }
            .branch { meta, sage_vcf ->

                def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.PAVE_VCF_NORMAL)

                runnable: Utils.hasTumorDnaBam(meta) && Utils.hasNormalDnaBam(meta) && sage_vcf && !has_existing
                skip: true
                    return meta
            }

        // Create process input channel
        // channel: [ meta_pave, sage_vcf ]
        ch_pave_germline_inputs = ch_sage_germline_inputs_sorted.runnable
            .map { meta, sage_vcf ->

                def meta_pave = [
                    key: meta.group_id,
                    id: meta.group_id,
                    sample_id: Utils.getTumorDnaSampleName(meta),
                ]

                return [meta_pave, sage_vcf]
            }

        // Run process
        GERMLINE(
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
            gnomad_resource,
        )

        ch_versions = ch_versions.mix(GERMLINE.out.versions)

        //
        // MODULE: PAVE somatic
        //
        // Select input sources and sort
        // channel: runnable: [ meta, sage_vcf ]
        // channel: skip: [ meta ]
        ch_sage_somatic_inputs_sorted = ch_sage_somatic_vcf
            .map { meta, sage_vcf, sage_tbi ->
                return [
                    meta,
                    Utils.selectCurrentOrExisting(sage_vcf, meta, Constants.INPUT.SAGE_VCF_TUMOR),
                ]
            }
            .branch { meta, sage_vcf ->

                def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.PAVE_VCF_TUMOR)

                runnable: Utils.hasTumorDnaBam(meta) && sage_vcf && !has_existing
                skip: true
                    return meta
            }

        // Create process input channel
        // channel: [ meta_pave, sage_vcf ]
        ch_pave_somatic_inputs = ch_sage_somatic_inputs_sorted.runnable
            .map { meta, sage_vcf ->

                def meta_pave = [
                    key: meta.group_id,
                    id: meta.group_id,
                    sample_id: Utils.getTumorDnaSampleName(meta),
                ]

                return [meta_pave, sage_vcf]
            }

        // Run process
        SOMATIC(
            ch_pave_somatic_inputs,
            genome_fasta,
            genome_version,
            genome_fai,
            sage_pon,
            pon_artefacts,
            segment_mappability,
            driver_gene_panel,
            ensembl_data_resources,
            gnomad_resource,
        )

        ch_versions = ch_versions.mix(SOMATIC.out.versions)

        // Set outputs, restoring original meta
        // channel: [ meta, gripss_vcf, gripss_tbi ]
        ch_somatic_out = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(SOMATIC.out.vcf, ch_inputs),
                ch_sage_somatic_inputs_sorted.skip.map { meta -> [meta, []] },
            )

        ch_germline_out = Channel.empty()
            .mix(
                WorkflowOncoanalyser.restoreMeta(GERMLINE.out.vcf, ch_inputs),
                ch_sage_germline_inputs_sorted.skip.map { meta -> [meta, []] },
            )

    emit:
        germline = ch_germline_out // channel: [ meta, pave_vcf ]
        somatic  = ch_somatic_out  // channel: [ meta, pave_vcf ]

        versions = ch_versions     // channel: [ versions.yml ]
}
