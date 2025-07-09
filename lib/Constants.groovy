class Constants {

    // NOTE(SW): the HMF reference data files are incompatible with hg19 due to different contig naming
    static List GENOMES_VERSION_37 = ['GRCh37_hmf', 'GRCh37']
    static List GENOMES_VERSION_38 = ['GRCh38_hmf', 'GRCh38', 'hg38']
    static List GENOMES_ALT        = ['GRCh38', 'hg38']

    static List GENOMES_SUPPORTED  = ['GRCh37_hmf', 'GRCh38_hmf']
    static List GENOMES_DEFINED    = Constants.GENOMES_VERSION_37 + Constants.GENOMES_VERSION_38

    static List PANELS_DEFINED     = ['tso500']


    static String HMF_DATA_37_PATH = 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.37_v2.1.0--1.tar.gz'
    static String HMF_DATA_38_PATH = 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.38_v2.1.0--1.tar.gz'

    static String TSO500_PANEL_37_PATH = 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/hmf_panel_resources.tso500.37_v2.0.0--3.tar.gz'
    static String TSO500_PANEL_38_PATH = 'https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/hmf_panel_resources.tso500.38_v2.0.0--3.tar.gz'

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
        CIDER,
        COBALT,
        CUPPA,
        ESVEE,
        ISOFOX,
        LILAC,
        LINX,
        NEO,
        ORANGE,
        PAVE,
        PEACH,
        PURPLE,
        REDUX,
        SAGE,
        SIGS,
        TEAL,
        VIRUSINTERPRETER,
    }

    static enum FileType {
        // Generic
        BAM,
        BAI,
        FASTQ,
        // Redux
        BAM_REDUX,
        REDUX_DUP_FREQ_TSV,
        REDUX_JITTER_TSV,
        REDUX_MS_TSV,
        // Process
        AMBER_DIR,
        BAMTOOLS_DIR,
        COBALT_DIR,
        ESVEE_VCF,
        ESVEE_VCF_TBI,
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
        LINX_PLOT_DIR,
        SAGE_DIR,
        PEACH_DIR,
    }

    static enum SampleType {
        TUMOR,
        NORMAL,
        TUMOR_NORMAL,
        DONOR,
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

        // Bams
        BAM_DNA_TUMOR: [
            FileType.BAM,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],

        BAM_DNA_NORMAL: [
            FileType.BAM,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],

        BAM_DNA_DONOR: [
            FileType.BAM,
            SampleType.DONOR,
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

        BAI_DNA_DONOR: [
            FileType.BAI,
            SampleType.DONOR,
            SequenceType.DNA,
        ],

        BAI_RNA_TUMOR: [
            FileType.BAI,
            SampleType.TUMOR,
            SequenceType.RNA,
        ],


        // REDUX
        BAM_REDUX_DNA_TUMOR: [
            FileType.BAM_REDUX,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],

        REDUX_DUP_FREQ_TSV_TUMOR: [
            FileType.REDUX_DUP_FREQ_TSV,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],

        REDUX_JITTER_TSV_TUMOR: [
            FileType.REDUX_JITTER_TSV,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],

        REDUX_MS_TSV_TUMOR: [
            FileType.REDUX_MS_TSV,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],

        BAM_REDUX_DNA_NORMAL: [
            FileType.BAM_REDUX,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],

        REDUX_DUP_FREQ_TSV_NORMAL: [
            FileType.REDUX_DUP_FREQ_TSV,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],

        REDUX_JITTER_TSV_NORMAL: [
            FileType.REDUX_JITTER_TSV,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],

        REDUX_MS_TSV_NORMAL: [
            FileType.REDUX_MS_TSV,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],

        BAM_REDUX_DNA_DONOR: [
            FileType.BAM_REDUX,
            SampleType.DONOR,
            SequenceType.DNA,
        ],

        REDUX_DUP_FREQ_TSV_DONOR: [
            FileType.REDUX_DUP_FREQ_TSV,
            SampleType.DONOR,
            SequenceType.DNA,
        ],

        REDUX_JITTER_TSV_DONOR: [
            FileType.REDUX_JITTER_TSV,
            SampleType.DONOR,
            SequenceType.DNA,
        ],

        REDUX_MS_TSV_DONOR: [
            FileType.REDUX_MS_TSV,
            SampleType.DONOR,
            SequenceType.DNA,
        ],


        // Other tools
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

        BAMTOOLS_DIR_TUMOR: [
            FileType.BAMTOOLS_DIR,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],
        BAMTOOLS_DIR_NORMAL: [
            FileType.BAMTOOLS_DIR,
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

        ESVEE_VCF_TUMOR: [
            FileType.ESVEE_VCF,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            SequenceType.DNA,
        ],
        ESVEE_VCF_TUMOR_TBI: [
            FileType.ESVEE_VCF_TBI,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            SequenceType.DNA,
        ],
        ESVEE_VCF_NORMAL: [
            FileType.ESVEE_VCF,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],
        ESVEE_VCF_NORMAL_TBI: [
            FileType.ESVEE_VCF_TBI,
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

        PEACH_DIR: [
            FileType.PEACH_DIR,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],

    ]
}
