class EmptyChannels {

    public static createEmptyChannel(ch, size) {
        // Input channel should have the structure: [meta, item_1, item_2, ...]
        return ch.map { d ->

            def meta = d[0]

            def output = [meta]
            size.times { output.add([]) }

            return output
        }
    }

    public static toolDir(ch) {
        return createEmptyChannel(ch, 1)
    }

    public static bamBai(ch) {
        return createEmptyChannel(ch, 2)
    }

    public static vcfTbi(ch) {
        return createEmptyChannel(ch, 2)
    }
}
