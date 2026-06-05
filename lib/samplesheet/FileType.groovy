package samplesheet

import util.Enums

public enum FileType {
    // Generic
    BAI,
    BAM,
    CRAI,
    CRAM,
    FASTQ,

    // REDUX
    BAM_REDUX,
    CRAM_REDUX,
    REDUX_TSV_DIR,

    // Other tools
    AMBER_DIR,
    BAMTOOLS_DIR,
    CHORD_DIR,
    COBALT_DIR,
    CUPPA_DIR,
    ESVEE_DIR,
    ISOFOX_DIR,
    LILAC_DIR,
    LINX_ANNO_DIR,
    LINX_PLOT_DIR,
    PAVE_DIR,
    PEACH_DIR,
    PURPLE_DIR,
    QSEE_DIR,
    SAGE_DIR,
    SAGE_APPEND_DIR,
    SAGE_VIS_DIR,
    SIGS_DIR,
    VIRUSINTERPRETER_DIR;

    public static FileType fromString(String string) {
        return Enums.getValidatedEnumFromString(string, FileType)
    }
}
