package samplesheet

import util.Enums

public enum InfoField {
    CANCER_TYPE,
    LANE,
    LIBRARY_ID,
    LONGITUDINAL_SAMPLE;

    public static InfoField fromString(String string) {
        return Enums.getValidatedEnumFromString(string, InfoField)
    }
}
