package pipeline

import util.Enums

public enum PurityEstimateMode {
    TARGETED,
    WGTS;

    public static PurityEstimateMode fromString(String string) {
        return Enums.getValidatedEnumFromString(string, PurityEstimateMode)
    }
}
