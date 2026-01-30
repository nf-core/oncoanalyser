class PlaceholderChannels {

    public static createPlaceholderChannel(ch, n_items) {
        // Input channel should have the structure `[meta, item_1, item_2, ...]` or `meta`
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
    public static final N_ITEMS_REDUX_TSVS = 4 // [ meta, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv ]

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
