//
// gene-utils region creation generates SAGE region files
//

import Constants
import Utils

include { GENE_UTILS_SAGE_REGIONS } from '../../../modules/local/gene_utils/sage_regions/main'

workflow GENE_UTILS_REGION_CREATION {
    take:
    // Input data
    driver_gene_panel      // channel: [mandatory] /path/to/driver_gene_panel

    // Reference data
    genome_version         // channel: [mandatory] genome version
    clinvar_annotations    // channel: [mandatory] /path/to/clinvar_annotations
    ensembl_data_resources // channel: [mandatory] /path/to/ensembl_data_resources/

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Run process
    GENE_UTILS_SAGE_REGIONS(
        driver_gene_panel,
        genome_version,
        clinvar_annotations,
        ensembl_data_resources,
    )

    ch_versions = ch_versions.mix(GENE_UTILS_SAGE_REGIONS.out.versions)

    emit:
    versions  = ch_versions // channel: [ versions.yml ]
}

