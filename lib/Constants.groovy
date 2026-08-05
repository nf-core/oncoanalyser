class Constants {

    // NOTE(SW): the HMF reference data files are incompatible with hg19 due to different contig naming
    static List GENOMES_VERSION_37 = ['GRCh37_hmf', 'GRCh37']
    static List GENOMES_VERSION_38 = ['GRCh38_hmf', 'GRCh38', 'hg38']
    static List GENOMES_ALT = ['GRCh38', 'hg38']

    static List GENOMES_SUPPORTED = ['GRCh37_hmf', 'GRCh38_hmf']
    static List GENOMES_DEFINED = Constants.GENOMES_VERSION_37 + Constants.GENOMES_VERSION_38

    static List PANELS_DEFINED = ['tso500']


    static String HMF_DATA_37_PATH = 'hartwig/pipeline_resources/hmf_pipeline_resources.37_v3.0.0--8.tar.gz'
    static String HMF_DATA_38_PATH = 'hartwig/pipeline_resources/hmf_pipeline_resources.38_v3.0.0--8.tar.gz'

    static String TSO500_PANEL_37_PATH = 'hartwig/panel_resources/hmf_panel_resources.tso500.37_v3.0.0--8.tar.gz'


    static Integer DEFAULT_ISOFOX_READ_LENGTH_WTS = 151
    static Integer DEFAULT_ISOFOX_READ_LENGTH_TARGETED = 93


    static enum RunMode {
        PANEL_RESOURCE_CREATION,
        PREPARE_REFERENCE,
        PURITY_ESTIMATE,
        TARGETED,
        WGTS,
    }

    static enum SequencingPlatform {
        ILLUMINA,
        SBX,
        ULTIMA,
    }

    static enum UmiType {
        KAPA,
        MSK,
        TSO500,
        TWIST,
    }

    static enum RefDataType {
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

    static List DEFAULT_EXCLUDED_PROCESSES = [] // For experimental tools

    static enum FileType {
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

    static enum SampleType {
        DONOR,
        NORMAL,
        TUMOR,
        TUMOR_NORMAL,
    }

    static enum SequenceType {
        DNA,
        DNA_RNA,
        RNA,
    }

    static enum InfoField {
        CANCER_TYPE,
        FLOWCELL,
        LANE,
        LIBRARY_ID,
        LONGITUDINAL_SAMPLE,
        GENERATE_REDUX_TSVS_ONLY,
        READ_GROUP_OVERRIDES,
    }

    static Map PLACEHOLDER_META = [meta_placeholder: null]

    static Map INPUT = [

        // Alignments
        // NOTE(SW): The ALN_DNA_* are only used in REDUX, where the ALN_REDUX can only run when they're set to have only TSVs generated
        ALN_DNA_TUMOR: [
            [FileType.ALN, FileType.ALN_REDUX],
            SampleType.TUMOR,
            SequenceType.DNA,
        ],

        ALN_DNA_NORMAL: [
            [FileType.ALN, FileType.ALN_REDUX],
            SampleType.NORMAL,
            SequenceType.DNA,
        ],

        ALN_DNA_DONOR: [
            [FileType.ALN, FileType.ALN_REDUX],
            SampleType.DONOR,
            SequenceType.DNA,
        ],

        ALN_RNA_TUMOR: [
            FileType.ALN,
            SampleType.TUMOR,
            SequenceType.RNA,
        ],

        IDX_DNA_TUMOR: [
            FileType.IDX,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],

        IDX_DNA_NORMAL: [
            FileType.IDX,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],

        IDX_DNA_DONOR: [
            FileType.IDX,
            SampleType.DONOR,
            SequenceType.DNA,
        ],

        IDX_RNA_TUMOR: [
            FileType.IDX,
            SampleType.TUMOR,
            SequenceType.RNA,
        ],


        // REDUX
        REDUX_DIR_TUMOR: [
            FileType.REDUX_DIR,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],

        REDUX_DIR_NORMAL: [
            FileType.REDUX_DIR,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],

        REDUX_DIR_DONOR: [
            FileType.REDUX_DIR,
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
        SAGE_APPEND_DIR_TUMOR: [
            FileType.SAGE_APPEND_DIR,
            SampleType.TUMOR,
            SequenceType.DNA_RNA,
        ],
        SAGE_APPEND_DIR_NORMAL: [
            FileType.SAGE_APPEND_DIR,
            SampleType.NORMAL,
            SequenceType.DNA_RNA,
        ],
        SAGE_PLOT_DIR_TUMOR: [
            FileType.SAGE_PLOT_DIR,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            SequenceType.DNA,
        ],

        PAVE_DIR_TUMOR: [
            FileType.PAVE_DIR,
            SampleType.TUMOR,
            SequenceType.DNA,
        ],
        PAVE_DIR_NORMAL: [
            FileType.PAVE_DIR,
            SampleType.NORMAL,
            SequenceType.DNA,
        ],

        ESVEE_DIR: [
            FileType.ESVEE_DIR,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            SequenceType.DNA,
        ],

        PURPLE_DIR: [
            FileType.PURPLE_DIR,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            SequenceType.DNA,
        ],

        QSEE_DIR: [
            FileType.QSEE_DIR,
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
