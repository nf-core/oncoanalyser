package refgenome

public enum SupportedGenome {

    // Genomes that are configured in igenomes.config or hmf_genomes.config
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

    public String defaultHmfDataPath() {
        return switch (this) {
            case GRCh37_hmf -> 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.37_v2.3.0--2.tar.gz'
            case GRCh38_hmf -> 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.38_v2.3.0--2.tar.gz'
        }
    }
}
