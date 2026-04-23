package pipeline

import refgenome.RefGenomeVersion

public enum SupportedPanel {
    TSO500,
    MSK,
    ONCOPANEL,
    PM_HAEM;

    public boolean hasConfiguredVersion(Map params, RefGenomeVersion version) {

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

}
