class PlaceholderChannels {

    public static createPlaceholderChannel(ch, size) {
        // Input channel should have the structure `[meta, item_1, item_2, ...]` or `meta`
        return ch.map { d ->
            def meta = (d instanceof List) ? d[0] : d

            def output = [meta]
            size.times { output.add([]) }

            return output
        }
    }

    public static toolDir(ch) {
        // Placeholder for: [ meta, dir ]
        return createPlaceholderChannel(ch, 1)
    }

    public static bamBai(ch) {
        // Placeholder for: [ meta, bam, bai ]
        return createPlaceholderChannel(ch, 2)
    }

    public static vcfTbi(ch) {
        // Placeholder for: [ meta, vcf, tbi ]
        return createPlaceholderChannel(ch, 2)
    }

    public static reduxTsvs(ch) {
        // Placeholder for: [ meta, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv ]
        return createPlaceholderChannel(ch, 4)
    }
}
