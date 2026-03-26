class RefGenome {

    public static enum Version {
        V37('37'),
        V38('38'),

        private final String numericName

        Version(String name) {
            this.numericName = name
        }

        public String getNumericName() { return this.numericName }

        public static List<String> getNumericNames() {
            return values().collect { version -> version.getNumericName() }
        }

        public static Version fromNumericName(String name) {
            switch (name) {
                case '37' -> V37
                case '38' -> V38
                default -> throw new IllegalArgumentException("Invalid ref genome version numeric name: ${name}")
            }
        }
    }

    public static enum Type {
        ALT,
        NO_ALT;

        public String getName() { return this.toString().toLowerCase() }

        public static List<String> getNames() {
            return values().collect { type -> type.toString() }
        }
    }

    public static enum SupportedGenome {

        // Genomes that are configured in igenomes.config or hmf_genomes.config
        GRCh37_hmf(Version.V37, Type.NO_ALT),
        GRCh38_hmf(Version.V38, Type.NO_ALT),

        private final Version version
        private final Type type

        SupportedGenome(Version version, Type type) {
            this.version = version
            this.type = type
        }

        public Version getVersion(){ return this.version }
        public Type getType(){ return this.type }

        public static SupportedGenome fromString(String string) {
            try {
                return valueOf(string)
            } catch (IllegalArgumentException e){
                return null
            }
        }

        public static List<String> getNames() {
            return values().collect{ variant -> variant.toString() }
        }
    }
}
