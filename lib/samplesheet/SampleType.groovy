package samplesheet

import util.Enums

public enum SampleType {
    DONOR,
    NORMAL,
    TUMOR,
    TUMOR_NORMAL;

    public static SampleType fromString(String string) {
        return Enums.getValidatedEnumFromString(string, SampleType)
    }
}
