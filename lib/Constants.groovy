class Constants {

    // NOTE(SW): the HMF reference data files are incompatible with hg19 due to different contig naming
    static List GENOMES_VERSION_37 = ['GRCh37_hmf', 'GRCh37']
    static List GENOMES_VERSION_38 = ['GRCh38_hmf', 'GRCh38', 'hg38']
    static List GENOMES_ALT        = ['GRCh38', 'hg38']

    static List GENOMES_SUPPORTED  = ['GRCh37_hmf', 'GRCh38_hmf']
    static List GENOMES_DEFINED    = Constants.GENOMES_VERSION_37 + Constants.GENOMES_VERSION_38

    static List PANELS_DEFINED     = ['hmf', 'tso500']


    static String HMF_DATA_37_PATH = 'https://pub-29f2e5b2b7384811bdbbcba44f8b5083.r2.dev/hmf_reference_data/hmftools/5.33_37--0.tar.gz'
    static String HMF_DATA_38_PATH = 'https://pub-29f2e5b2b7384811bdbbcba44f8b5083.r2.dev/hmf_reference_data/hmftools/5.33_38--0.tar.gz'


    static String HMF_PANEL_38_PATH = 'https://pub-29f2e5b2b7384811bdbbcba44f8b5083.r2.dev/hmf_reference_data/panels/hmf_5.33_38--0.tar.gz'

    static String TSO500_PANEL_37_PATH = 'https://pub-29f2e5b2b7384811bdbbcba44f8b5083.r2.dev/hmf_reference_data/panels/tso500_5.33_37--0.tar.gz'
    static String TSO500_PANEL_38_PATH = 'https://pub-29f2e5b2b7384811bdbbcba44f8b5083.r2.dev/hmf_reference_data/panels/tso500_5.33_38--0.tar.gz'


    static String VIRUSBREAKENDDB_PATH = 'https://pub-29f2e5b2b7384811bdbbcba44f8b5083.r2.dev/virusbreakend/virusbreakenddb_20210401.tar.gz'

    static String HLA_SLICE_BED_GRCH38_ALT_PATH = 'https://pub-29f2e5b2b7384811bdbbcba44f8b5083.r2.dev/umccr_reference_data/other/hla_slice/grch38_alt.plus_homologous.bed'

    static enum RunType {
        TUMOR_NORMAL,
        TUMOR_ONLY,
    }

    static enum RunMode {
        PANEL,
        WGS,
        WGTS,
        WTS,
    }

    static enum Process {
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
        // Process
        AMBER_DIR,
        BAMTOOLS,
        COBALT_DIR,
        GRIDSS_VCF,
        GRIPSS_VCF,
        GRIPSS_UNFILTERED_VCF,
        ISOFOX_DIR,
        LILAC_DIR,
        LINX_ANNO_DIR,
        PAVE_VCF,
        PURPLE_DIR,
        SAGE_VCF,
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
        TARGETTED,
        WGS,
        WGTS,
        WTS,
    }

    static Map PLACEHOLDER_META = [meta_placeholder: null]
    static List PLACEHOLDER_OPTIONAL_CHANNEL = []

    static Map INPUT = [

        ISOFOX_DIR: [
            FileType.ISOFOX_DIR,
            SampleType.TUMOR,
            SequenceType.WTS
        ],

        AMBER_DIR: [
            FileType.AMBER_DIR,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            [SequenceType.TARGETTED, SequenceType.WGS]
        ],
        COBALT_DIR: [
            FileType.COBALT_DIR,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            [SequenceType.TARGETTED, SequenceType.WGS]
        ],

        BAMTOOLS_TUMOR: [
            FileType.BAMTOOLS,
            SampleType.TUMOR,
            [SequenceType.TARGETTED, SequenceType.WGS],
        ],
        BAMTOOLS_NORMAL: [
            FileType.BAMTOOLS,
            SampleType.NORMAL,
            SequenceType.WGS,
        ],

        FLAGSTAT_TUMOR: [
            FileType.FLAGSTAT,
            SampleType.TUMOR,
            [SequenceType.TARGETTED, SequenceType.WGS],
        ],
        FLAGSTAT_NORMAL: [
            FileType.FLAGSTAT,
            SampleType.NORMAL,
            SequenceType.WGS
        ],

        SAGE_VCF_TUMOR: [
            FileType.SAGE_VCF,
            SampleType.TUMOR,
            [SequenceType.TARGETTED, SequenceType.WGS],
        ],
        SAGE_VCF_NORMAL: [
            FileType.SAGE_VCF,
            SampleType.NORMAL,
            SequenceType.WGS
        ],
        SAGE_DIR_TUMOR: [
            FileType.SAGE_DIR,
            SampleType.TUMOR,
            [SequenceType.TARGETTED, SequenceType.WGS],
        ],
        SAGE_DIR_NORMAL: [
            FileType.SAGE_DIR,
            SampleType.NORMAL,
            SequenceType.WGS
        ],

        PAVE_VCF_TUMOR: [
            FileType.PAVE_VCF,
            SampleType.TUMOR,
            [SequenceType.TARGETTED, SequenceType.WGS],
        ],
        PAVE_VCF_NORMAL: [
            FileType.PAVE_VCF,
            SampleType.NORMAL,
            SequenceType.WGS
        ],

        GRIDSS_VCF: [
            FileType.GRIDSS_VCF,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            [SequenceType.TARGETTED, SequenceType.WGS],
        ],

        GRIPSS_VCF_TUMOR: [
            FileType.GRIPSS_VCF,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            [SequenceType.TARGETTED, SequenceType.WGS],
        ],
        GRIPSS_VCF_NORMAL: [
            FileType.GRIPSS_VCF,
            SampleType.NORMAL,
            SequenceType.WGS
        ],
        GRIPSS_UNFILTERED_VCF_TUMOR: [
            FileType.GRIPSS_UNFILTERED_VCF,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            [SequenceType.TARGETTED, SequenceType.WGS],
        ],
        GRIPSS_UNFILTERED_VCF_NORMAL: [
            FileType.GRIPSS_UNFILTERED_VCF,
            SampleType.NORMAL,
            SequenceType.WGS
        ],

        PURPLE_DIR: [
            FileType.PURPLE_DIR,
            [SampleType.TUMOR, SampleType.TUMOR_NORMAL],
            [SequenceType.TARGETTED, SequenceType.WGS],
        ],

        LINX_PLOT_DIR_TUMOR: [
            FileType.LINX_PLOT_DIR,
            SampleType.TUMOR,
            [SequenceType.TARGETTED, SequenceType.WGS],
        ],
        LINX_ANNO_DIR_TUMOR: [
            FileType.LINX_ANNO_DIR,
            SampleType.TUMOR,
            [SequenceType.TARGETTED, SequenceType.WGS],
        ],
        LINX_ANNO_DIR_NORMAL: [
            FileType.LINX_ANNO_DIR,
            SampleType.NORMAL,
            SequenceType.WGS
        ],

        CHORD_DIR: [
            FileType.CHORD_DIR,
            SampleType.TUMOR,
            SequenceType.WGS
        ],
        SIGS_DIR: [
            FileType.SIGS_DIR,
            SampleType.TUMOR,
            SequenceType.WGS
        ],
        LILAC_DIR: [
            FileType.LILAC_DIR,
            [SampleType.TUMOR, SampleType.NORMAL, SampleType.TUMOR_NORMAL],
            [SequenceType.TARGETTED, SequenceType.WGS, SequenceType.WGTS],
        ],

        VIRUSINTERPRETER_DIR: [
            FileType.VIRUSINTERPRETER_DIR,
            SampleType.TUMOR,
            SequenceType.WGS
        ],

        CUPPA_DIR: [
            FileType.CUPPA_DIR,
            SampleType.TUMOR,
            [SequenceType.WGS, SequenceType.WTS, SequenceType.WGTS],
        ],

    ]
}
