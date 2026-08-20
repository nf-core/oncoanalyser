//
// Shared enumerations and derived type-level data for the oncoanalyser pipeline.
// These live here (not lib/*.groovy) so they can be explicitly included, which avoids
// a name collision between the `Process` enum and java.lang.Process.
//

enum RunMode {
    PANEL_RESOURCE_CREATION,
    PREPARE_REFERENCE,
    PURITY_ESTIMATE,
    TARGETED,
    WGTS,
}

enum SequencingPlatform {
    ILLUMINA,
    SBX,
    ULTIMA,
}

enum UmiType {
    KAPA,
    MSK,
    TSO500,
    TWIST,
}

enum RefDataType {
    // Compound types
    TARGETED,
    WGS,
    WTS,

    // Individual types
    BWAMEM2_INDEX,
    DICT,
    DNA_ALIGNMENT,
    FAI,
    FASTA,
    GRIDSS_INDEX,
    HMFTOOLS,
    IMG,
    PANEL,
    RNA_ALIGNMENT,
    STAR_INDEX,
}

enum Process {
    ALIGNMENT,
    AMBER,
    BAMTOOLS,
    CHORD,
    CIDER,
    COBALT,
    CUPPA,
    ESVEE,
    ISOFOX,
    LILAC,
    LINX,
    MULTIQC,
    NEO,
    ORANGE,
    PAVE,
    PEACH,
    PURPLE,
    QSEE,
    REDUX,
    SAGE,
    SAGE_APPEND,
    SAGE_VISUALISER,
    SIGS,
    TEAL,
    VIRUSINTERPRETER,
    WISP,
}

enum FileType {
    // Input only
    FASTQ,
    BAM,
    CRAM,
    BAM_REDUX,
    CRAM_REDUX,
    BAI,
    CRAI,

    // Generic (internal representation)
    ALN,
    ALN_REDUX,
    IDX,

    // Process
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
    REDUX_DIR,
    SAGE_APPEND_DIR,
    SAGE_DIR,
    SAGE_PLOT_DIR,
    SIGS_DIR,
    VIRUSINTERPRETER_DIR,
}

enum SampleType {
    DONOR,
    LONGITUDINAL,
    NORMAL,
    TUMOR,
    TUMOR_NORMAL,
}

enum SequenceType {
    DNA,
    DNA_RNA,
    RNA,
}

enum InfoField {
    CANCER_TYPE,
    FLOWCELL,
    LANE,
    LIBRARY_ID,
    LONGITUDINAL_SAMPLE,
    GENERATE_REDUX_TSVS_ONLY,
    READ_GROUP_OVERRIDES,
}

// Directories that are produced/consumed once per case rather than per sample. These live on
// CaseRecord.directories instead of a SampleRecord.files entry.
def getCaseLevelDirs() {
    return [
        FileType.AMBER_DIR,
        FileType.COBALT_DIR,
        FileType.ESVEE_DIR,
        FileType.PURPLE_DIR,
        FileType.QSEE_DIR,
        FileType.SAGE_PLOT_DIR,
        FileType.LINX_PLOT_DIR,
        FileType.LILAC_DIR,
        FileType.CHORD_DIR,
        FileType.SIGS_DIR,
        FileType.VIRUSINTERPRETER_DIR,
        FileType.CUPPA_DIR,
        FileType.PEACH_DIR,
    ] as Set
}
