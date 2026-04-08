//
// PAVE PON creation prepares the panel-specific small variant artefact resource
//

include { PAVE_PON_PANEL_CREATION } from '../../../modules/local/pave/pon_creation/main'


workflow PAVE_PON_CREATION {
    take:
    // Sample data
    ch_sage_dir_somatic // channel: [mandatory] [ meta, sage_dir ]

    // Reference data
    genome_version      // channel: [mandatory] genome version

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Create process input channel
    // channel: [ [sage_vcf, ...], [sage_tbi, ...] ]
    ch_pave_inputs = ch_sage_dir_somatic
        .map { meta, sage_dir ->

            def (sage_vcf, sage_tbi) = Inputs.resolveSageVcfWithTbi(sage_dir, meta, SampleSheetFields.SampleType.TUMOR)

            return [ sage_vcf, sage_tbi ]
        }
        .collect(flat: false)
        .map { d -> d.transpose() }

    // Run process
    PAVE_PON_PANEL_CREATION(
        ch_pave_inputs,
        genome_version,
    )

    ch_versions = ch_versions.mix(PAVE_PON_PANEL_CREATION.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
