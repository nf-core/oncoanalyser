package refgenome

public enum RefGenomeVersion {
    V37('37'),
    V38('38'),

    private final String numericName

    RefGenomeVersion(String name) {
        this.numericName = name
    }

    public String getNumericName() { return this.numericName }

    public static List<String> getNumericNames() {
        return values().collect { version -> version.getNumericName() }
    }

    public static RefGenomeVersion fromNumericName(String name) {
        return switch (name) {
            case '37' -> V37
            case '38' -> V38
            default -> throw new IllegalArgumentException("Invalid ref genome version numeric name: ${name}")
        }
    }

    public SupportedGenome getSupportedGenome(){
        return switch (this) {
            case V37 -> SupportedGenome.GRCh37_hmf
            case V38 -> SupportedGenome.GRCh38_hmf
        }
    }
}
