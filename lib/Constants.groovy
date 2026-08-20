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
        LONGITUDINAL,
        NORMAL,
        TUMOR,
        TUMOR_NORMAL,
    }

    // Directories that are produced/consumed once per case rather than per sample. These live on
    // CaseRecord.directories instead of a SampleRecord.files entry.
    static Set CASE_LEVEL_DIRS = [
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
    ]

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

}
