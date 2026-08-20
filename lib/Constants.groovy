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

    static List DEFAULT_EXCLUDED_PROCESSES = [] // For experimental tools

    static Map PLACEHOLDER_META = [meta_placeholder: null]

}
