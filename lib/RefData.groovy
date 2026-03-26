class RefData {

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

    public static String getDefaultHmfDataPath(RefGenome.Version version) {

        return switch (version) {
            case RefGenome.Version.V37 -> 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.37_v2.3.0--2.tar.gz'
            case RefGenome.Version.V38 -> 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.38_v2.3.0--2.tar.gz'
            default -> throw new IllegalArgumentException("Invalid ref genome version: ${version}")
        }
    }

    public static enum SupportedPanel {
        tso500,
        msk,
        oncopanel,
        pm_haem;

        public boolean hasConfiguredVersion(Map params, RefGenome.Version version) {

            def panel_versions = params.panel_data_paths[this.toString()]
            def versioned_panel_data_paths = panel_versions[version.getNumericName()]

            return versioned_panel_data_paths != null
        }

        public static SupportedPanel fromString(String string) {
            try {
                return valueOf(string)
            } catch (IllegalArgumentException e){
                return null
            }
        }

        public static List<String> getNames() {
            return values().collect{ panel -> panel.toString() }
        }

    }

    public static String getDefaultPanelDataPath(SupportedPanel panel, RefGenome.Version version){

        def path_map = [
            [(SupportedPanel.tso500), (RefGenome.Version.V37)]: 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/hmf_panel_resources.tso500.37_v2.3.0--2.tar.gz',
            [(SupportedPanel.tso500), (RefGenome.Version.V38)]: 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/hmf_panel_resources.tso500.38_v2.3.0--2.tar.gz',
        ]

        def panel_key = [panel, version]

        if (!path_map.containsKey(panel_key))
            throw new NoSuchElementException("No default panel data path defined for panel(${panel}) refGenomeVersion(${version.getNumericName()})")

        return path_map[panel_key]
    }
}
