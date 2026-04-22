package channels

class PlaceholderChannels {

    private static Object createPlaceholderChannel(Object ch, int size) {
        // Given a channel containing metadata with the structure `[meta, ...]` or atomic: `meta`,
        // Create a channel with the form e.g. `[meta, [], []]` if size=2
        return ch.map { d ->
            def meta = (d instanceof List) ? d[0] : d

            def output = [meta]
            size.times { output.add([]) }

            return output
        }
    }

    public static final int N_ITEMS_TOOL_DIR = 1 // [ meta, dir ]
    public static final int N_ITEMS_BAM_BAI = 2 // [ meta, bam, bai ]
    public static final int N_ITEMS_VCF_TBI = 2 // [ meta, vcf, tbi ]
    public static final int N_ITEMS_REDUX_TSV_LIST = 1 // [ meta, [redux_tsv, ...] ]

    public static Object toolDir(Object ch) {
        return createPlaceholderChannel(ch, N_ITEMS_TOOL_DIR)
    }

    public static Object bamBai(Object ch) {
        return createPlaceholderChannel(ch, N_ITEMS_BAM_BAI)
    }

    public static Object vcfTbi(Object ch) {
        return createPlaceholderChannel(ch, N_ITEMS_VCF_TBI)
    }

    public static Object reduxTsvs(Object ch) {
        return createPlaceholderChannel(ch, N_ITEMS_REDUX_TSV_LIST)
    }
}
