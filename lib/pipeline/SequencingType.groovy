package pipeline

import util.Enums

public enum SequencingType {
    ILLUMINA,
    SBX,
    ULTIMA;

    public static SequencingType fromString(String string) {
        return Enums.getValidatedEnumFromString(string, SequencingType, false)
    }
}
