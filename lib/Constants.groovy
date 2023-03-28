class Constants {

    // NOTE(SW): the HMF reference data files are incompatible with hg19 due to different contig naming
    static List GENOMES_VERSION_37       = ['GRCh37_hmf', 'GRCh37']
    static List GENOMES_VERSION_38       = ['GRCh38_hmf', 'GRCh38', 'hg38']
    static List GENOMES_ALT              = ['GRCh38', 'hg38']

    static String HMF_DATA_37_PATH = 'https://pub-29f2e5b2b7384811bdbbcba44f8b5083.r2.dev/hmf_reference_data/repacks/5.31_37_0.0.1.tar.gz'
    static String HMF_DATA_38_PATH = 'https://pub-29f2e5b2b7384811bdbbcba44f8b5083.r2.dev/hmf_reference_data/repacks/5.32+dev1_38_0.0.1.tar.gz'

    static String VIRUSBREAKENDDB_PATH = 'https://pub-29f2e5b2b7384811bdbbcba44f8b5083.r2.dev/virusbreakend/virusbreakenddb_20210401.tar.gz'

    static enum PipelineMode {
        FULL,
        MANUAL,
        GRIDSS_PURPLE_LINX,
        CUPPA,
    }

    static enum Process {
        AMBER,
        BAMTOOLS,
        CHORD,
        COBALT,
        CUPPA,
        GRIDSS,
        GRIPSS,
        ISOFOX,
        LILAC,
        LINX,
        ORANGE,
        PAVE,
        PEACH,
        PROTECT,
        PURPLE,
        SAGE,
        SIGS,
        SVPREP,
        VIRUSINTERPRETER,
    }

    static enum FileType {
        // Generic
        BAM,
        // Process
        AMBER_DIR,
        BAMTOOLS_TXT,
        COBALT_DIR,
        GRIDSS_VCF,
        GRIPSS_HARD_VCF,
        GRIPSS_SOFT_VCF,
        ISOFOX_DIR,
        LILAC_DIR,
        LINX_ANNO_DIR,
        PAVE_VCF,
        PURPLE_DIR,
        SAGE_VCF,
        VIRUSINTERPRETER_TSV,
        // ORANGE specific
        CHORD_PREDICTION,
        CUPPA_CSV,
        CUPPA_FEATURE_PLOT,
        CUPPA_SUMMARY_PLOT,
        FLAGSTAT,
        LINX_PLOT_DIR,
        PEACH_TSV,
        PROTECT_TSV,
        SAGE_BQR,
        SAGE_COVERAGE,
    }

    static enum SampleType {
        TUMOR,
        NORMAL,
        TUMOR_NORMAL,
    }

    static enum SequenceType {
        WGS,
        WGTS,
        WTS,
    }

    static Map META_PLACEHOLDER = [meta_placeholder: null]

    static Map INPUT = [
        ISOFOX_DIR:             [FileType.ISOFOX_DIR,           SampleType.TUMOR,        SequenceType.WTS],

        AMBER_DIR:              [FileType.AMBER_DIR,            SampleType.TUMOR_NORMAL, SequenceType.WGS],
        COBALT_DIR:             [FileType.COBALT_DIR,           SampleType.TUMOR_NORMAL, SequenceType.WGS],

        BAMTOOLS_TXT_TUMOR:     [FileType.BAMTOOLS_TXT,         SampleType.TUMOR,        SequenceType.WGS],
        BAMTOOLS_TXT_NORMAL:    [FileType.BAMTOOLS_TXT,         SampleType.NORMAL,       SequenceType.WGS],

        SAGE_VCF_TUMOR:         [FileType.SAGE_VCF,             SampleType.TUMOR,        SequenceType.WGS],
        SAGE_VCF_NORMAL:        [FileType.SAGE_VCF,             SampleType.NORMAL,       SequenceType.WGS],
        SAGE_BQR_TUMOR:         [FileType.SAGE_BQR,             SampleType.TUMOR,        SequenceType.WGS],
        SAGE_BQR_NORMAL:        [FileType.SAGE_BQR,             SampleType.NORMAL,       SequenceType.WGS],
        SAGE_COVERAGE:          [FileType.SAGE_COVERAGE,        SampleType.NORMAL,       SequenceType.WGS],

        PAVE_VCF_TUMOR:         [FileType.PAVE_VCF,             SampleType.TUMOR,        SequenceType.WGS],
        PAVE_VCF_NORMAL:        [FileType.PAVE_VCF,             SampleType.NORMAL,       SequenceType.WGS],

        GRIDSS_VCF:             [FileType.GRIDSS_VCF,           SampleType.TUMOR_NORMAL, SequenceType.WGS],

        GRIPSS_HARD_VCF_TUMOR:  [FileType.GRIPSS_HARD_VCF,      SampleType.TUMOR,        SequenceType.WGS],
        GRIPSS_SOFT_VCF_TUMOR:  [FileType.GRIPSS_SOFT_VCF,      SampleType.TUMOR,        SequenceType.WGS],
        GRIPSS_HARD_VCF_NORMAL: [FileType.GRIPSS_HARD_VCF,      SampleType.NORMAL,       SequenceType.WGS],
        GRIPSS_SOFT_VCF_NORMAL: [FileType.GRIPSS_SOFT_VCF,      SampleType.NORMAL,       SequenceType.WGS],

        PURPLE_DIR:             [FileType.PURPLE_DIR,           SampleType.TUMOR_NORMAL, SequenceType.WGS],

        LINX_PLOT_DIR_TUMOR:    [FileType.LINX_PLOT_DIR,        SampleType.TUMOR,        SequenceType.WGS],
        LINX_ANNO_DIR_TUMOR:    [FileType.LINX_ANNO_DIR,        SampleType.TUMOR,        SequenceType.WGS],
        LINX_ANNO_DIR_NORMAL:   [FileType.LINX_ANNO_DIR,        SampleType.NORMAL,       SequenceType.WGS],

        CHORD_PREDICTION:       [FileType.CHORD_PREDICTION,     SampleType.TUMOR,        SequenceType.WGS],
        LILAC_DIR:              [FileType.LILAC_DIR,            SampleType.NORMAL,       SequenceType.WGS],
        PEACH_TSV:              [FileType.PEACH_TSV,            SampleType.NORMAL,       SequenceType.WGS],
        PROTECT_TSV:            [FileType.PROTECT_TSV,          SampleType.TUMOR,        SequenceType.WGS],

        VIRUSINTERPRETER_TSV:   [FileType.VIRUSINTERPRETER_TSV, SampleType.TUMOR,        SequenceType.WGS],

        CUPPA_CSV:              [FileType.CUPPA_CSV,            SampleType.TUMOR,        SequenceType.WGTS],
        CUPPA_FEATURE_PLOT:     [FileType.CUPPA_FEATURE_PLOT,   SampleType.TUMOR,        SequenceType.WGTS],
        CUPPA_SUMMARY_PLOT:     [FileType.CUPPA_SUMMARY_PLOT,   SampleType.TUMOR,        SequenceType.WGTS],
    ]
}
