package refgenome

public enum SupportedGenome {

    // Genomes that are configured in igenomes.config or hmf_genomes.config
    GRCh37_hmf(RefGenomeVersion.V37, RefGenomeType.NO_ALT),
    GRCh38_hmf(RefGenomeVersion.V38, RefGenomeType.NO_ALT),

    private final RefGenomeVersion version
    private final RefGenomeType type

    SupportedGenome(RefGenomeVersion version, RefGenomeType type) {
        this.version = version
        this.type = type
    }

    public RefGenomeVersion getVersion(){ return this.version }
    public RefGenomeType getType(){ return this.type }

    public static SupportedGenome fromString(String string) {
        try {
            return valueOf(string)
        } catch (IllegalArgumentException e){
            return null
        }
    }
}
