package refgenome

public enum SupportedGenome {

    GRCh37_hmf(RefGenomeVersion.V37, RefGenomeType.NO_ALT),
    GRCh38_hmf(RefGenomeVersion.V38, RefGenomeType.NO_ALT);

    private final RefGenomeVersion genomeVersion
    private final RefGenomeType genomeType

    SupportedGenome(RefGenomeVersion genomeVersion, RefGenomeType genomeType) {
        this.genomeVersion = genomeVersion
        this.genomeType = genomeType
    }

    public RefGenomeVersion genomeVersion(){ return this.genomeVersion }
    public RefGenomeType genomeType(){ return this.genomeType }

    public static SupportedGenome fromString(String string) {
        try {
            return valueOf(string)
        } catch (IllegalArgumentException e){
            return null
        }
    }
}
