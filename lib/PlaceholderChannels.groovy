class PlaceholderChannels {

    private static createPlaceholderChannel(ch, n_items) {
        // Given a channel containing metadata with the structure `[meta, ...]` or atomic: `meta`,
        // Create a channel with the form e.g. `[meta, [], []]` if n_items=2
        return ch.map { d ->
            def meta = (d instanceof List) ? d[0] : d

            def output = [meta]
            n_items.times { output.add([]) }

            return output
        }
    }

    public static final N_ITEMS_TOOL_DIR = 1   // [ meta, dir ]
    public static final N_ITEMS_BAM_BAI = 2    // [ meta, bam, bai ]
    public static final N_ITEMS_VCF_TBI = 2    // [ meta, vcf, tbi ]
    public static final N_ITEMS_REDUX_TSVS = 3 // [ meta, bqr_tsv, jitter_tsv, ms_tsv ]

    public static toolDir(ch) {
        return createPlaceholderChannel(ch, N_ITEMS_TOOL_DIR)
    }

    public static bamBai(ch) {
        return createPlaceholderChannel(ch, N_ITEMS_BAM_BAI)
    }

    public static vcfTbi(ch) {
        return createPlaceholderChannel(ch, N_ITEMS_VCF_TBI)
    }

    public static reduxTsvs(ch) {
        return createPlaceholderChannel(ch, N_ITEMS_REDUX_TSVS)
    }
}
