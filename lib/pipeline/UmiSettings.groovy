package pipeline

public enum UmiSettings {

    TSO500(
        new FastpArgs(false, NO_UMI_LOCATION, NO_UMI_LENGTH, NO_UMI_SKIP),
        new FastqToolsArgs(false, NO_DUPLEX_UMI_DELIM),
        new ReduxArgs(true, "+")
    ),

    TWIST(
        new FastpArgs(true, "per_read", 7, 0),
        new FastqToolsArgs(false, NO_DUPLEX_UMI_DELIM),
        new ReduxArgs(true, "_")
    ),

    MSK(
        new FastpArgs(false, NO_UMI_LOCATION, NO_UMI_LENGTH, NO_UMI_SKIP),
        new FastqToolsArgs(true, "+"),
        new ReduxArgs(true, NO_DUPLEX_UMI_DELIM)
    ),

    NONE(
        new FastpArgs(false, NO_UMI_LOCATION, NO_UMI_LENGTH, NO_UMI_SKIP),
        new FastqToolsArgs(false, NO_DUPLEX_UMI_DELIM),
        new ReduxArgs(false, NO_DUPLEX_UMI_DELIM)
    );

    private final FastpArgs fastpArgs
    private final FastqToolsArgs fastqToolsArgs
    private final ReduxArgs reduxArgs

    UmiSettings(FastpArgs fastpArgs, FastqToolsArgs fastqToolsArgs, ReduxArgs reduxArgs) {
        this.fastpArgs = fastpArgs
        this.fastqToolsArgs = fastqToolsArgs
        this.reduxArgs = reduxArgs
    }

    public FastpArgs fastpArgs(){ return this.fastpArgs }
    public FastqToolsArgs fastqToolsArgs(){ return this.fastqToolsArgs }
    public ReduxArgs reduxArgs(){ return this.reduxArgs }

    public static UmiSettings fromSupportedPanel(SupportedPanel supportedPanel) {
        return switch (supportedPanel) {
            case SupportedPanel.TSO500 -> TSO500
            case SupportedPanel.MSK -> MSK
            case SupportedPanel.PM_HAEM -> TWIST
            default -> NONE
        }
    }

    public static record FastpArgs(boolean requireUmiStripping, String umiLocation, int umiLength, int umiSkip) {}
    public static record FastqToolsArgs(boolean requireUmiStripping, String umiDelim) {}
    public static record ReduxArgs(boolean enableUmiProcessing, String duplexUmiDelim) {}

    private static final String NO_UMI_LOCATION = ""
    private static final int NO_UMI_LENGTH = 0
    private static final int NO_UMI_SKIP = -1
    private static final String NO_DUPLEX_UMI_DELIM = ""
}