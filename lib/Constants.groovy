class Constants {

    // NOTE(SW): the HMF reference data files are incompatible with hg19 due to different contig naming
    static List GENOMES_VERSION_37 = ['GRCh37_hmf', 'GRCh37']
    static List GENOMES_VERSION_38 = ['GRCh38_hmf', 'GRCh38', 'hg38']
    static List GENOMES_ALT        = ['GRCh38', 'hg38']

    static List GENOMES_SUPPORTED  = ['GRCh37_hmf', 'GRCh38_hmf']
    static List GENOMES_DEFINED    = Constants.GENOMES_VERSION_37 + Constants.GENOMES_VERSION_38

    static List PANELS_DEFINED     = ['tso500']


    static String HMF_DATA_37_PATH = 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/5.34_37--2.tar.gz'
    static String HMF_DATA_38_PATH = 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/5.34_38--2.tar.gz'


    static String TSO500_PANEL_37_PATH = 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/tso500_5.34_37--1.tar.gz'
    static String TSO500_PANEL_38_PATH = 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/tso500_5.34_38--1.tar.gz'


    static String VIRUSBREAKENDDB_PATH = 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/virusbreakend/virusbreakenddb_20210401.tar.gz'

    static String HLA_SLICE_BED_GRCH38_ALT_PATH = 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/other/hla_slice/grch38_alt.plus_homologous.bed'


    static Integer DEFAULT_ISOFOX_READ_LENGTH_WTS = 151
    static Integer DEFAULT_ISOFOX_READ_LENGTH_TARGETED = 93


    static enum RunMode {
        TARGETED,
        WGTS,
    }

    static enum Process {
        ALIGNMENT,
        AMBER,
        BAMTOOLS,
        CHORD,
        COBALT,
        CUPPA,
        FLAGSTAT,
        GRIDSS,
        GRIPSS,
        ISOFOX,
        LILAC,
        LINX,
        MARKDUPS,
        ORANGE,
        PAVE,
        PURPLE,
        SAGE,
        SIGS,
        VIRUSINTERPRETER,
    }

    static enum FileType {
        // Generic
        BAM,
        BAM_MARKDUPS,
        BAI,
        FASTQ,
        // Process
        AMBER_DIR,
        BAMTOOLS,
        COBALT_DIR,
        GRIDSS_VCF,
        GRIDSS_VCF_TBI,
        GRIPSS_VCF,
        GRIPSS_VCF_TBI,
        GRIPSS_UNFILTERED_VCF,
        GRIPSS_UNFILTERED_VCF_TBI,
        ISOFOX_DIR,
        LILAC_DIR,
        LINX_ANNO_DIR,
        PAVE_VCF,
        PURPLE_DIR,
        SAGE_VCF,
        SAGE_VCF_TBI,
        SAGE_APPEND_VCF,
        VIRUSINTERPRETER_DIR,
        // ORANGE specific
        CHORD_DIR,
        SIGS_DIR,
        CUPPA_DIR,
        FLAGSTAT,
        LINX_PLOT_DIR,
        SAGE_DIR,
    }

    static enum SampleType {
        TUMOR,
        NORMAL,
        TUMOR_NORMAL,
    }

    static enum SequenceType {
        DNA,
        RNA,
        DNA_RNA,
    }

    static enum InfoField {
        CANCER_TYPE,
        LANE,
        LIBRARY_ID,
    }

    static Map PLACEHOLDER_META = [meta_placeholder: null]
    static List PLACEHOLDER_OPTIONAL_CHANNEL = []

    static Map INPUT = [

        BAM_DNA_TUMOR: [
            FileType.BAM,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],

        BAM_MARKDUPS_DNA_TUMOR: [
            FileType.BAM_MARKDUPS,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],

        BAM_DNA_NORMAL: [
            FileType.BAM,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],

        BAM_MARKDUPS_DNA_NORMAL: [
            FileType.BAM_MARKDUPS,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],

        BAM_RNA_TUMOR: [
            FileType.BAM,
            SampleType.TUMOR,
            SequenceType.RNA,
        ],

        BAI_DNA_TUMOR: [
            FileType.BAI,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],

        BAI_DNA_NORMAL: [
            FileType.BAI,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],

        BAI_RNA_TUMOR: [
            FileType.BAI,
            SampleType.TUMOR,
            SequenceType.RNA,
        ],

        ISOFOX_DIR: [
            FileType.ISOFOX_DIR,
            SampleType.TUMOR,
            SequenceType.RNA,
        ],

        AMBER_DIR: [
            FileType.AMBER_DIR,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            SequenceType.DNA,
        ],
        COBALT_DIR: [
            FileType.COBALT_DIR,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            SequenceType.DNA,
        ],

        BAMTOOLS_TUMOR: [
            FileType.BAMTOOLS,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],
        BAMTOOLS_NORMAL: [
            FileType.BAMTOOLS,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],

        FLAGSTAT_TUMOR: [
            FileType.FLAGSTAT,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],
        FLAGSTAT_NORMAL: [
            FileType.FLAGSTAT,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],

        SAGE_VCF_TUMOR: [
            FileType.SAGE_VCF,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],
        SAGE_VCF_NORMAL: [
            FileType.SAGE_VCF,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],
        SAGE_VCF_TBI_TUMOR: [
            FileType.SAGE_VCF_TBI,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],
        SAGE_VCF_TBI_NORMAL: [
            FileType.SAGE_VCF_TBI,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],
        SAGE_DIR_TUMOR: [
            FileType.SAGE_DIR,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],
        SAGE_DIR_NORMAL: [
            FileType.SAGE_DIR,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],
        SAGE_APPEND_VCF_TUMOR: [
            FileType.SAGE_APPEND_VCF,
            SampleType.TUMOR,
            SequenceType.DNA_RNA,
        ],
        SAGE_APPEND_VCF_NORMAL: [
            FileType.SAGE_APPEND_VCF,
            SampleType.NORMAL,
            SequenceType.DNA_RNA,
        ],

        PAVE_VCF_TUMOR: [
            FileType.PAVE_VCF,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],
        PAVE_VCF_NORMAL: [
            FileType.PAVE_VCF,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],

        GRIDSS_VCF: [
            FileType.GRIDSS_VCF,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            SequenceType.DNA,
        ],

        GRIPSS_VCF_TUMOR: [
            FileType.GRIPSS_VCF,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            SequenceType.DNA,
        ],
        GRIPSS_VCF_TUMOR_TBI: [
            FileType.GRIPSS_VCF_TBI,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            SequenceType.DNA,
        ],
        GRIPSS_VCF_NORMAL: [
            FileType.GRIPSS_VCF,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],
        GRIPSS_VCF_NORMAL_TBI: [
            FileType.GRIPSS_VCF_TBI,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],
        GRIPSS_UNFILTERED_VCF_TUMOR: [
            FileType.GRIPSS_UNFILTERED_VCF,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            SequenceType.DNA,
        ],
        GRIPSS_UNFILTERED_VCF_TUMOR_TBI: [
            FileType.GRIPSS_UNFILTERED_VCF_TBI,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            SequenceType.DNA,
        ],
        GRIPSS_UNFILTERED_VCF_NORMAL: [
            FileType.GRIPSS_UNFILTERED_VCF,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],

        PURPLE_DIR: [
            FileType.PURPLE_DIR,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            SequenceType.DNA,
        ],

        LINX_PLOT_DIR_TUMOR: [
            FileType.LINX_PLOT_DIR,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],
        LINX_ANNO_DIR_TUMOR: [
            FileType.LINX_ANNO_DIR,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],
        LINX_ANNO_DIR_NORMAL: [
            FileType.LINX_ANNO_DIR,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],

        CHORD_DIR: [
            FileType.CHORD_DIR,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],
        SIGS_DIR: [
            FileType.SIGS_DIR,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],
        LILAC_DIR: [
            FileType.LILAC_DIR,
            [SampleType.TUMOR, SampleType.NORMAL, SampleType.TUMOR_NORMAL],
            [SequenceType.DNA, SequenceType.DNA_RNA],
        ],

        VIRUSINTERPRETER_DIR: [
            FileType.VIRUSINTERPRETER_DIR,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],

        CUPPA_DIR: [
            FileType.CUPPA_DIR,
            SampleType.TUMOR,
            [SequenceType.DNA, SequenceType.RNA, SequenceType.DNA_RNA],
        ],

    ]
}
