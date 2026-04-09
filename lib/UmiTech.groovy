public enum UmiTech {

    TSO500(
        new FastpParams(false, "per_read", 15, -1),
        new ReduxParams(true, "+")
    ),

    TWIST(
        new FastpParams(true, "per_read", 7, 0),
        new ReduxParams(true, "_")
    ),

    NONE(
        new FastpParams(false, "", 0, -1),
        new ReduxParams(false, "")
    );

    private final FastpParams fastpParams
    private final ReduxParams reduxParams

    UmiTech (FastpParams fastpParams, ReduxParams reduxParams) {
        this.fastpParams = fastpParams
        this.reduxParams = reduxParams
    }

    public FastpParams fastpParams(){ return this.fastpParams }
    public ReduxParams reduxParams(){ return this.reduxParams }

    public static UmiTech fromSupportedPanel(RefData.SupportedPanel supportedPanel) {
        return switch (supportedPanel) {
            case RefData.SupportedPanel.tso500 -> TSO500
            case RefData.SupportedPanel.pm_haem -> TWIST
            default -> NONE
        }
    }

    public static record FastpParams(boolean requireUmiStripping, String umiLocation, int umiLength, int umiSkip) {}

    public static record ReduxParams(boolean enableUmiProcessing, String duplexUmiDelim) {}
}