package panel

import refgenome.RefGenomeVersion

public enum SupportedPanel {

    TSO500_GRCH37(
        Panel.TSO500,
        RefGenomeVersion.V37,
        'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/hmf_panel_resources.tso500.37_v2.3.0--2.tar.gz'
    ),

    MSK_GRCH37(
        Panel.MSK,
        RefGenomeVersion.V37,
        null
    );

    private final Panel panel
    private final RefGenomeVersion genomeVersion
    private final String defaultDataPath

    SupportedPanel(Panel panel, RefGenomeVersion genomeVersion, String defaultDataPath) {
        this.panel = panel
        this.genomeVersion = genomeVersion
        this.defaultDataPath = defaultDataPath
    }

    public Panel panel(){ return this.panel }
    public RefGenomeVersion refGenomeVersion() { return this.genomeVersion }

    public String defaultDataPath() {

        if(this.defaultDataPath == null)
            throw new NoSuchElementException("No default data path defined for panel(${panel}) refGenomeVersion($genomeVersion)")

        return this.defaultDataPath
    }

    public static SupportedPanel from(Panel panel, RefGenomeVersion genomeVersion) {

        def matchedSupportedPanel = values().find { supportedPanel ->
            panel == supportedPanel.panel && genomeVersion == supportedPanel.genomeVersion
        }

        if(!matchedSupportedPanel) {
            return null
        }

        return matchedSupportedPanel
    }

    public static SupportedPanel from(String panelStr, String genomeVersionStr) {
        def panel = Panel.fromString(panelStr)
        def genomeVersion = RefGenomeVersion.fromNumericName(genomeVersionStr)
        return from(panel, genomeVersion)
    }

    public static List<String> listAll(){
        return values().collect { supportedPanel ->
            return "panel(${supportedPanel.panel}) refGenomeVersion(${supportedPanel.genomeVersion.getNumericName()})"
        }
    }
}
