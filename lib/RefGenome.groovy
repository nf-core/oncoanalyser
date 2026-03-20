class RefGenome {

    // NOTE(SW): the HMF reference data files are incompatible with hg19 due to different contig naming
    public static final List<String> GENOMES_VERSION_37 = ['GRCh37_hmf', 'GRCh37']
    public static final List<String> GENOMES_VERSION_38 = ['GRCh38_hmf', 'GRCh38', 'hg38']
    public static final List<String> GENOMES_ALT = ['GRCh38', 'hg38']

    public static final List<String> GENOMES_SUPPORTED = ['GRCh37_hmf', 'GRCh38_hmf']
    public static final List<String> GENOMES_DEFINED = GENOMES_VERSION_37 + GENOMES_VERSION_38

    public static enum Version {
        V37('37'),
        V38('38');

        private final String name

        Version(String name) {
            this.name = name
        }

        public String getName() { return this.name }
    }

    public static enum Type {
        ALT,
        NO_ALT;

        public String getName() { return this.toString().toLowerCase() }
    }
}
