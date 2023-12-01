include { MARKDUPS } from '../../modules/local/markdups/main'

workflow READ_PROCESSING {
    take:
        // Sample data
        ch_inputs   // channel: [mandatory] [ meta ]
        ch_dna_bams // channel: [mandatory] [ meta, bam_dna ]
        ch_rna_bams // channel: [mandatory] [ meta, bam_rna ]

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // NOTE(SW): channel operations will be required to configure MarkDups for individual samples

        // MarkDups
        // TODO(SW): implement inputs
        ch_markdups_inputs = Channel.of([[id: 'foo'], []])
        MARKDUPS(
            ch_markdups_inputs,
            // TODO(SW): configuration
            // TODO(SW): reference files
        )
        // TODO(SW): implement outputs
        ch_markduplicates_dna_out = Channel.empty()
        ch_markduplicates_rna_out = Channel.empty()

    emit:
        dna       = ch_markduplicates_dna_out // channel: [ meta, bam_dna ]
        rna       = ch_markduplicates_rna_out // channel: [ meta, bam_rna ]

        versions  = ch_versions               // channel: [ versions.yml ]
}
