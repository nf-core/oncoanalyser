class Constants {

    // NOTE(SW): the HMF reference data files are incompatible with hg19 due to different contig naming
    static List GENOMES_VERSION_37       = ['GRCh37_hmf', 'GRCh37']
    static List GENOMES_VERSION_38       = ['GRCh38_hmf', 'GRCh38', 'hg38']
    static List GENOMES_ALT              = ['GRCh38', 'hg38']

    static String HMF_DATA_37_BASE = 's3://umccr-research-dev/stephen/hmf_reference_data/2_repacked/hmf_data_bundle_5.29_37_0.0.1/'
    static String HMF_DATA_38_BASE = 's3://umccr-research-dev/stephen/hmf_reference_data/2_repacked/hmf_data_bundle_5.29_38_0.0.1/'

    static String VIRUSBREAKENDDB_PATH = 's3://virusbreakend/virusbreakenddb_20210401.tar.gz'

    static enum PipelineMode {
        FULL,
        MANUAL,
        GRIDSS_PURPLE_LINX,
        CUPPA,
    }

    static enum Process {
        AMBER,
        CHORD,
        COBALT,
        COLLECTWGSMETRICS,
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
        VIRUSINTERPRETER,
    }

    static enum DataType {
        TUMOR,
        NORMAL,
        TUMOR_NORMAL,
    }

    static enum FileType {
        // Generic
        BAM_WGS,
        BAM_WTS,
        // Process
        AMBER_DIR,
        COBALT_DIR,
        COLLECTWGSMETRICS,
        GRIDSS_VCF,
        GRIPSS_HARD_VCF,
        GRIPSS_SOFT_VCF,
        ISOFOX_DIR,
        LINX_ANNO_DIR,
        LILAC_DIR,
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
}
