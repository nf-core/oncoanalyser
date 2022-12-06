class Constants {

    // NOTE(SW): the HMF reference data files are incompatible with hg19 due to different contig naming
    static List genomes_version_37       = ['GRCh37_hmf', 'GRCh37']
    static List genomes_version_38       = ['GRCh38', 'hg38']
    static List genomes_version_38_noalt = ['GRCh38_hmf']

    static String hmf_reference_data_37_bundle_path = 'PLACEHOLDER_hmf_reference_data_37_bundle_path'
    static String hmf_reference_data_38_bundle_path = 'PLACEHOLDER_hmf_reference_data_38_bundle_path'

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
        PAVE,
        PURPLE,
        SAGE,
        SIGS,
        SVPREP,
        TEAL,
        VIRUSINTERPRETER,
    }

    static public java.util.Map HMF_DATA_PATHS = [
        AMBER_LOCI:                   'amber/GermlineHetPon.vcf.gz',
        COBALT_GC_PROFILE:            'cobalt/DiploidRegions.bed.gz',
        CUPPA:                        'cuppa/',
        SV_PREP_BLACKLIST:            'svprep/sv_prep_blacklist.bed',
        GRIDSS_BLACKLIST:             'gridss/ENCFF356LFX.bed.gz',
        GRIDSS_BREAKEND_PON:          'gridss/gridss_pon_single_breakend.bed.gz',
        GRIDSS_BREAKPOINT_PON:        'gridss/gridss_pon_breakpoint.bedpe.gz',
        REPEAT_MASKER_FILE:           'gridss/repeat_masker.fa.out.gz',
        ISOFOX_EXP_COUNTS:            'isofox/read_151_exp_counts.csv',
        ISOFOX_EXP_GC_RATIOS:         'isofox/read_100_exp_gc_ratios.csv',
        LINX_FRAGILE_SITES:           'linx/fragile_sites_hmf.csv',
        LINX_LINES:                   'linx/line_elements.csv',
        SAGE_BLACKLIST_BED:           'sage/KnownBlacklist.germline.bed',
        SAGE_BLACKLIST_VCF:           'sage/KnownHotspots.somatic.vcf.gz',
        SAGE_CODING_PANEL:            'sage/ActionableCodingPanel.bed.gz',
        SAGE_HIGH_CONFIDENCE:         'sage/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz',
        SAGE_KNOWN_HOTSPOTS_GERMLINE: 'sage/KnownBlacklist.germline.vcf.gz',
        SAGE_KNOWN_HOTSPOTS_SOMATIC:  'sage/KnownHotspots.somatic.vcf.gz',
        SAGE_PON_FILE:                'sage/SageGermlinePon.98x.tsv.gz',
        CLINVAR_VCF:                  'sage/clinvar.vcf.gz',
        SIGS_SIGNATURES:              'sigs/snv_cosmic_signatures.csv',
        LILAC_RESOURCE_DIR:           'lilac/',
        VIRUSBREAKENDDB:              'virusbreakend/',
        VIRUSINTERPRETER_TAXONOMY:    'virusinterpreter/taxonomy_db.tsv',
        VIRUSINTERPRETER_REPORTING:   'virusinterpreter/virus_reporting_db.tsv',
        PURPLE_GERMLINE_DEL:          'purple/cohort_germline_del_freq.csv',
        DRIVER_GENE_PANEL:            'gene_panel/DriverGenePanel.tsv',
        ENSEMBL_DATA_DIR:             'ensembl_data_cache/',
        KNOWN_FUSION_DATA:            'known_fusions/known_fusion_data.csv',
        KNOWN_FUSIONS:                'known_fusions/known_fusions.bedpe',
        MAPPABILITY_BED:              'mappability/mappability_150.bed.gz',
    ]
}
