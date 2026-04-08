class SampleMeta {

    public static enum FileType {
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
        SIGS_DIR,
        VIRUSINTERPRETER_DIR;

        public static FileType fromString(String string) {
            return Enums.getValidatedEnumFromString(string, FileType)
        }
    }

    public static enum SampleType {
        DONOR,
        NORMAL,
        TUMOR,
        TUMOR_NORMAL;

        public static SampleType fromString(String string) {
            return Enums.getValidatedEnumFromString(string, SampleType)
        }
    }

    public static enum SequenceType {
        DNA,
        DNA_RNA,
        RNA;

        public static SequenceType fromString(String string) {
            return Enums.getValidatedEnumFromString(string, SequenceType)
        }
    }

    public static enum InfoField {
        CANCER_TYPE,
        LANE,
        LIBRARY_ID,
        LONGITUDINAL_SAMPLE;

        public static InfoField fromString(String string) {
            return Enums.getValidatedEnumFromString(string, InfoField)
        }
    }

    public static final Map INPUT = [

        // Bams
        BAM_DNA_TUMOR:  [FileType.BAM, SampleType.TUMOR, SequenceType.DNA],
        BAM_DNA_NORMAL: [FileType.BAM, SampleType.NORMAL, SequenceType.DNA],
        BAM_DNA_DONOR:  [FileType.BAM, SampleType.DONOR, SequenceType.DNA],
        BAM_RNA_TUMOR:  [FileType.BAM, SampleType.TUMOR, SequenceType.RNA],

        BAI_DNA_TUMOR:  [FileType.BAI, SampleType.TUMOR, SequenceType.DNA],
        BAI_DNA_NORMAL: [FileType.BAI, SampleType.NORMAL, SequenceType.DNA],
        BAI_DNA_DONOR:  [FileType.BAI, SampleType.DONOR, SequenceType.DNA],
        BAI_RNA_TUMOR:  [FileType.BAI, SampleType.TUMOR, SequenceType.RNA],

        // REDUX
        BAM_REDUX_DNA_TUMOR:  [FileType.BAM_REDUX, SampleType.TUMOR, SequenceType.DNA],
        BAM_REDUX_DNA_NORMAL: [FileType.BAM_REDUX, SampleType.NORMAL, SequenceType.DNA],
        BAM_REDUX_DNA_DONOR:  [FileType.BAM_REDUX, SampleType.DONOR, SequenceType.DNA],

        REDUX_TSV_DIR_TUMOR:  [FileType.REDUX_TSV_DIR, SampleType.TUMOR, SequenceType.DNA],
        REDUX_TSV_DIR_NORMAL: [FileType.REDUX_TSV_DIR, SampleType.NORMAL, SequenceType.DNA],
        REDUX_TSV_DIR_DONOR:  [FileType.REDUX_TSV_DIR, SampleType.DONOR, SequenceType.DNA],

        // Other tools
        AMBER_DIR: [FileType.AMBER_DIR, [SampleType.TUMOR, SampleType.TUMOR_NORMAL], SequenceType.DNA],

        BAMTOOLS_DIR_TUMOR: [FileType.BAMTOOLS_DIR, SampleType.TUMOR, SequenceType.DNA],
        BAMTOOLS_DIR_NORMAL: [FileType.BAMTOOLS_DIR, SampleType.NORMAL, SequenceType.DNA],

        CHORD_DIR: [FileType.CHORD_DIR, SampleType.TUMOR, SequenceType.DNA],

        COBALT_DIR: [FileType.COBALT_DIR, [SampleType.TUMOR, SampleType.TUMOR_NORMAL], SequenceType.DNA],

        CUPPA_DIR: [FileType.CUPPA_DIR, SampleType.TUMOR, [SequenceType.DNA, SequenceType.RNA, SequenceType.DNA_RNA]],

        ESVEE_DIR: [FileType.ESVEE_DIR, [SampleType.TUMOR, SampleType.TUMOR_NORMAL], SequenceType.DNA],

        ISOFOX_DIR: [FileType.ISOFOX_DIR, SampleType.TUMOR, SequenceType.RNA],

        LILAC_DIR: [
            FileType.LILAC_DIR,
            [SampleType.TUMOR, SampleType.NORMAL, SampleType.TUMOR_NORMAL],
            [SequenceType.DNA, SequenceType.DNA_RNA],
        ],

        LINX_PLOT_DIR_TUMOR:  [FileType.LINX_PLOT_DIR, SampleType.TUMOR, SequenceType.DNA],
        LINX_ANNO_DIR_TUMOR:  [FileType.LINX_ANNO_DIR, SampleType.TUMOR, SequenceType.DNA],
        LINX_ANNO_DIR_NORMAL: [FileType.LINX_ANNO_DIR, SampleType.NORMAL, SequenceType.DNA],

        SAGE_APPEND_DIR_TUMOR:  [FileType.SAGE_APPEND_DIR, SampleType.TUMOR, SequenceType.DNA_RNA],
        SAGE_APPEND_DIR_NORMAL: [FileType.SAGE_APPEND_DIR, SampleType.NORMAL, SequenceType.DNA_RNA],

        SAGE_DIR_TUMOR:  [FileType.SAGE_DIR, SampleType.TUMOR, SequenceType.DNA],
        SAGE_DIR_NORMAL: [FileType.SAGE_DIR, SampleType.NORMAL, SequenceType.DNA],

        PAVE_DIR_TUMOR:  [FileType.PAVE_DIR, SampleType.TUMOR, SequenceType.DNA],
        PAVE_DIR_NORMAL: [FileType.PAVE_DIR, SampleType.NORMAL, SequenceType.DNA],

        PEACH_DIR: [FileType.PEACH_DIR, SampleType.NORMAL, SequenceType.DNA],

        PURPLE_DIR: [FileType.PURPLE_DIR, [SampleType.TUMOR, SampleType.TUMOR_NORMAL], SequenceType.DNA],

        QSEE_DIR: [FileType.QSEE_DIR, [SampleType.TUMOR, SampleType.TUMOR_NORMAL], SequenceType.DNA],

        SIGS_DIR: [FileType.SIGS_DIR, SampleType.TUMOR, SequenceType.DNA],

        VIRUSINTERPRETER_DIR: [FileType.VIRUSINTERPRETER_DIR, SampleType.TUMOR, SequenceType.DNA],
    ]
}
