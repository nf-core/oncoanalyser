package samplesheet

import util.Enums

public enum SequenceType {
    DNA,
    DNA_RNA,
    RNA;

    public static SequenceType fromString(String string) {
        return Enums.getValidatedEnumFromString(string, SequenceType)
    }
}
