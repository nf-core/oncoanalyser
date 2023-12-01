include { BWA_MEM2 } from '../../modules/local/bwa/mem2/main'
include { STAR     } from '../../modules/local/star/main'

workflow READ_ALIGNMENT {
    take:
        // Sample data
        ch_inputs // channel: [mandatory] [ meta ]

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // STAR
        // TODO(SW): implement inputs
        ch_star_inputs = Channel.of([[id: 'foo'], []])
        STAR(
          ch_star_inputs,
          // TODO(SW): include reference files
        )
        // TODO(SW): implement outputs
        ch_star_outputs = Channel.empty()

        // BWA MEM2
        // TODO(SW): implement inputs
        ch_bwa_inputs = Channel.of([[id: 'foo'], []])
        BWA_MEM2(
          ch_bwa_inputs,
          // TODO(SW): include reference files
        )
        // TODO(SW): implement outputs
        ch_bwa_outputs = Channel.empty()

    emit:
        dna       = ch_bwa_outputs  // channel: [ meta, bam_dna ]
        rna       = ch_star_outputs // channel: [ meta, bam_rna ]

        versions  = ch_versions     // channel: [ versions.yml ]
}
