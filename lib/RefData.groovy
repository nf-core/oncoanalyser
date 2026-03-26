class RefData {

    public static final String HMF_DATA_37_PATH = 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.37_v2.3.0--2.tar.gz'
    public static final String HMF_DATA_38_PATH = 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.38_v2.3.0--2.tar.gz'

    public static final String TSO500_PANEL_37_PATH = 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/hmf_panel_resources.tso500.37_v2.3.0--2.tar.gz'
    public static final String TSO500_PANEL_38_PATH = 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/hmf_panel_resources.tso500.38_v2.3.0--2.tar.gz'

    public static final List<String> PANELS_DEFINED = ['tso500', 'msk', 'oncopanel', 'pm_haem']

    public static enum Type {
        // Compound types
        TARGETED,
        WGS,
        WTS,

        // Individual types
        BWAMEM2_INDEX,
        DICT,
        DNA_ALIGNMENT,
        FAI,
        FASTA,
        GRIDSS_INDEX,
        HMFTOOLS,
        IMG,
        PANEL,
        RNA_ALIGNMENT,
        STAR_INDEX,
    }
}
