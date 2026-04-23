package pipeline

import util.Enums

public enum UmiType {

    TSO500(
        fastpUmiDisabled(),
        fastqToolsUmiDisabled(),
        new ReduxArgs(true, "+")
    ),

    TWIST(
        new FastpArgs(true, "per_read", 7, 0),
        fastqToolsUmiDisabled(),
        new ReduxArgs(true, "_")
    ),

    MSK(
        fastpUmiDisabled(),
        new FastqToolsArgs(true, "+"),
        new ReduxArgs(true, NO_DUPLEX_UMI_DELIM)
    ),

    NONE(
        fastpUmiDisabled(),
        fastqToolsUmiDisabled(),
        reduxUmiDisabled()
    );

    private final FastpArgs fastpArgs
    private final FastqToolsArgs fastqToolsArgs
    private final ReduxArgs reduxArgs

    UmiType(FastpArgs fastpArgs, FastqToolsArgs fastqToolsArgs, ReduxArgs reduxArgs) {
        this.fastpArgs = fastpArgs
        this.fastqToolsArgs = fastqToolsArgs
        this.reduxArgs = reduxArgs
    }

    public static UmiType fromString(String type){
        return Enums.getValidatedEnumFromString(type, UmiType)
    }

    public static UmiType fromSupportedPanel(SupportedPanel supportedPanel) {
        return switch (supportedPanel) {
            case SupportedPanel.TSO500 -> TSO500
            case SupportedPanel.MSK -> MSK
            case SupportedPanel.PM_HAEM -> TWIST
            default -> NONE
        }
    }

    public FastpArgs fastpArgs(){ return this.fastpArgs }
    public FastqToolsArgs fastqToolsArgs(){ return this.fastqToolsArgs }
    public ReduxArgs reduxArgs(){ return this.reduxArgs }

    public static record FastpArgs(boolean umiProcessingEnabled, String umiLocation, int umiLength, int umiSkip) {}
    public static record FastqToolsArgs(boolean umiProcessingEnabled, String umiDelim) {}
    public static record ReduxArgs(boolean umiProcessingEnabled, String duplexUmiDelim) {}

    private static final String NO_UMI_LOCATION = ""
    private static final int NO_UMI_LENGTH = 0
    private static final int NO_UMI_SKIP = -1
    private static final String NO_UMI_DELIM = ""
    private static final String NO_DUPLEX_UMI_DELIM = ""

    private static FastpArgs fastpUmiDisabled() { return new FastpArgs(false, NO_UMI_LOCATION, NO_UMI_LENGTH, NO_UMI_SKIP) }
    private static FastqToolsArgs fastqToolsUmiDisabled() { return new FastqToolsArgs(false, NO_UMI_DELIM) }
    private static ReduxArgs reduxUmiDisabled() { return new ReduxArgs(false, NO_DUPLEX_UMI_DELIM) }

}