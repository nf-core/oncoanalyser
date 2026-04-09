public enum UmiSettings {

    TSO500(
        new FastpArgs(true, "per_read", 15, 0),
        new ReduxArgs(true, "+")
    ),

    TWIST(
        new FastpArgs(true, "per_read", 7, 0),
        new ReduxArgs(true, "_")
    ),

    NONE(
        new FastpArgs(false, "", 0, -1),
        new ReduxArgs(false, "")
    );

    private final FastpArgs fastpArgs
    private final ReduxArgs reduxArgs

    UmiSettings(FastpArgs fastpArgs, ReduxArgs reduxArgs) {
        this.fastpArgs = fastpArgs
        this.reduxArgs = reduxArgs
    }

    public FastpArgs fastpArgs(){ return this.fastpArgs }
    public ReduxArgs reduxArgs(){ return this.reduxArgs }

    public static UmiSettings fromSupportedPanel(RefData.SupportedPanel supportedPanel) {
        return switch (supportedPanel) {
            case RefData.SupportedPanel.TSO500 -> TSO500
            case RefData.SupportedPanel.PM_HAEM -> TWIST
            default -> NONE
        }
    }

    public static record FastpArgs(boolean requireUmiStripping, String umiLocation, int umiLength, int umiSkip) {}

    public static record ReduxArgs(boolean enableUmiProcessing, String duplexUmiDelim) {}
}