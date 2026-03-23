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
        REDUX_BQR_TSV,
        REDUX_DUP_FREQ_TSV,
        REDUX_JITTER_TSV,
        REDUX_MS_TABLE_TSV,
        REDUX_MSI_TSV,

        // Process
        AMBER_DIR,
        BAMTOOLS_DIR,
        COBALT_DIR,
        ESVEE_DIR,
        ISOFOX_DIR,
        LILAC_DIR,
        LINX_ANNO_DIR,
        PAVE_DIR,
        PURPLE_DIR,
        SAGE_VCF,
        SAGE_VCF_TBI,
        SAGE_APPEND_DIR,
        VIRUSINTERPRETER_DIR,

        // ORANGE specific
        CHORD_DIR,
        CUPPA_DIR,
        LINX_PLOT_DIR,
        PEACH_DIR,
        SAGE_DIR,
        SIGS_DIR,
    }

    public static enum SampleType {
        DONOR,
        NORMAL,
        TUMOR,
        TUMOR_NORMAL,
    }

    public static enum SequenceType {
        DNA,
        DNA_RNA,
        RNA,
    }

    public static enum InfoField {
        CANCER_TYPE,
        LANE,
        LIBRARY_ID,
        LONGITUDINAL_SAMPLE,
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

        REDUX_BQR_TSV_TUMOR:  [FileType.REDUX_BQR_TSV, SampleType.TUMOR, SequenceType.DNA],
        REDUX_BQR_TSV_NORMAL: [FileType.REDUX_BQR_TSV, SampleType.NORMAL, SequenceType.DNA],
        REDUX_BQR_TSV_DONOR:  [FileType.REDUX_BQR_TSV, SampleType.DONOR, SequenceType.DNA],

        REDUX_DUP_FREQ_TSV_TUMOR:  [FileType.REDUX_DUP_FREQ_TSV, SampleType.TUMOR, SequenceType.DNA],
        REDUX_DUP_FREQ_TSV_NORMAL: [FileType.REDUX_DUP_FREQ_TSV, SampleType.NORMAL, SequenceType.DNA],
        REDUX_DUP_FREQ_TSV_DONOR:  [FileType.REDUX_DUP_FREQ_TSV, SampleType.DONOR, SequenceType.DNA],

        REDUX_JITTER_TSV_TUMOR:  [FileType.REDUX_JITTER_TSV, SampleType.TUMOR, SequenceType.DNA],
        REDUX_JITTER_TSV_NORMAL: [FileType.REDUX_JITTER_TSV, SampleType.NORMAL, SequenceType.DNA],
        REDUX_JITTER_TSV_DONOR:  [FileType.REDUX_JITTER_TSV, SampleType.DONOR, SequenceType.DNA],

        REDUX_MS_TABLE_TSV_TUMOR:  [FileType.REDUX_MS_TABLE_TSV, SampleType.TUMOR, SequenceType.DNA],
        REDUX_MS_TABLE_TSV_NORMAL: [FileType.REDUX_MS_TABLE_TSV, SampleType.NORMAL, SequenceType.DNA],
        REDUX_MS_TABLE_TSV_DONOR:  [FileType.REDUX_MS_TABLE_TSV, SampleType.DONOR, SequenceType.DNA],

        REDUX_MSI_TSV_TUMOR:  [FileType.REDUX_MSI_TSV, SampleType.TUMOR, SequenceType.DNA],
        REDUX_MSI_TSV_NORMAL: [FileType.REDUX_MSI_TSV, SampleType.NORMAL, SequenceType.DNA],
        REDUX_MSI_TSV_DONOR:  [FileType.REDUX_MSI_TSV, SampleType.DONOR, SequenceType.DNA],

        // Other tools
        ISOFOX_DIR: [FileType.ISOFOX_DIR, SampleType.TUMOR, SequenceType.RNA],

        AMBER_DIR: [FileType.AMBER_DIR, [SampleType.TUMOR, SampleType.TUMOR_NORMAL], SequenceType.DNA],

        COBALT_DIR: [FileType.COBALT_DIR, [SampleType.TUMOR, SampleType.TUMOR_NORMAL], SequenceType.DNA],

        BAMTOOLS_DIR_TUMOR: [FileType.BAMTOOLS_DIR, SampleType.TUMOR, SequenceType.DNA],
        BAMTOOLS_DIR_NORMAL: [FileType.BAMTOOLS_DIR, SampleType.NORMAL, SequenceType.DNA],

        SAGE_VCF_TUMOR:      [FileType.SAGE_VCF, SampleType.TUMOR, SequenceType.DNA],
        SAGE_VCF_NORMAL:     [FileType.SAGE_VCF, SampleType.NORMAL, SequenceType.DNA],
        SAGE_VCF_TBI_TUMOR:  [FileType.SAGE_VCF_TBI, SampleType.TUMOR, SequenceType.DNA],
        SAGE_VCF_TBI_NORMAL: [FileType.SAGE_VCF_TBI, SampleType.NORMAL, SequenceType.DNA],
        SAGE_DIR_TUMOR:      [FileType.SAGE_DIR, SampleType.TUMOR, SequenceType.DNA],
        SAGE_DIR_NORMAL:     [FileType.SAGE_DIR, SampleType.NORMAL, SequenceType.DNA],

        SAGE_APPEND_DIR_TUMOR:  [FileType.SAGE_APPEND_DIR, SampleType.TUMOR, SequenceType.DNA_RNA],
        SAGE_APPEND_DIR_NORMAL: [FileType.SAGE_APPEND_DIR, SampleType.NORMAL, SequenceType.DNA_RNA],

        PAVE_DIR_TUMOR:  [FileType.PAVE_DIR, SampleType.TUMOR, SequenceType.DNA],
        PAVE_DIR_NORMAL: [FileType.PAVE_DIR, SampleType.NORMAL, SequenceType.DNA],

        ESVEE_DIR: [FileType.ESVEE_DIR, [SampleType.TUMOR, SampleType.TUMOR_NORMAL], SequenceType.DNA],

        PURPLE_DIR: [FileType.PURPLE_DIR, [SampleType.TUMOR, SampleType.TUMOR_NORMAL], SequenceType.DNA],

        LINX_PLOT_DIR_TUMOR:  [FileType.LINX_PLOT_DIR, SampleType.TUMOR, SequenceType.DNA],
        LINX_ANNO_DIR_TUMOR:  [FileType.LINX_ANNO_DIR, SampleType.TUMOR, SequenceType.DNA],
        LINX_ANNO_DIR_NORMAL: [FileType.LINX_ANNO_DIR, SampleType.NORMAL, SequenceType.DNA],

        CHORD_DIR: [FileType.CHORD_DIR, SampleType.TUMOR, SequenceType.DNA],

        SIGS_DIR: [FileType.SIGS_DIR, SampleType.TUMOR, SequenceType.DNA],

        LILAC_DIR: [
            FileType.LILAC_DIR,
            [SampleType.TUMOR, SampleType.NORMAL, SampleType.TUMOR_NORMAL],
            [SequenceType.DNA, SequenceType.DNA_RNA],
        ],

        VIRUSINTERPRETER_DIR: [FileType.VIRUSINTERPRETER_DIR, SampleType.TUMOR, SequenceType.DNA],

        CUPPA_DIR: [FileType.CUPPA_DIR, SampleType.TUMOR, [SequenceType.DNA, SequenceType.RNA, SequenceType.DNA_RNA]],

        PEACH_DIR: [FileType.PEACH_DIR, SampleType.NORMAL, SequenceType.DNA],
    ]
}
