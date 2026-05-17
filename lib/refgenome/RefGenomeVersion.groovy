package refgenome

public enum RefGenomeVersion {
    V37('37', 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.37_v2.3.0--2.tar.gz'),
    V38('38', 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.38_v2.3.0--2.tar.gz'),

    private final String numericName
    private final String defaultHmfDataPath

    RefGenomeVersion(String name, String defaultHmfDataPath) {
        this.numericName = name
        this.defaultHmfDataPath = defaultHmfDataPath
    }

    public String getNumericName() { return this.numericName }
    public String defaultHmfDataPath() { return this.defaultHmfDataPath }

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
}
