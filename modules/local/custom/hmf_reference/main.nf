process HMF_REFERENCE {
    container 'public.ecr.aws/ubuntu/ubuntu:20.04_stable'

    input:
    path hmf_bundle

    output:
    path 'output'      , emit: dir
    path 'versions.yml', emit: versions

    script:
    """
    # if tar.gz
    # tar -xzvf ${hmf_bundle} --strip-components 1 -C output/
    # else
    # ln -s ${hmf_bundle} output/

    sleep 5

    touch versions.yml

    mkdir -p output/amber/
    touch output/amber/GermlineHetPon.vcf.gz

    mkdir -p output/cobalt/
    touch output/cobalt/DiploidRegions.bed.gz

    mkdir -p output/cuppa/

    mkdir -p output/svprep/
    touch output/svprep/sv_prep_blacklist.bed

    mkdir -p output/gridss/
    touch output/gridss/ENCFF356LFX.bed.gz
    touch output/gridss/gridss_pon_single_breakend.bed.gz
    touch output/gridss/gridss_pon_breakpoint.bedpe.gz
    touch output/gridss/repeat_masker.fa.out.gz

    mkdir -p output/isofox/
    touch output/isofox/read_151_exp_counts.csv
    touch output/isofox/read_100_exp_gc_ratios.csv

    mkdir -p output/linx/
    touch output/linx/fragile_sites_hmf.csv
    touch output/linx/line_elements.csv

    mkdir -p output/sage/
    touch output/sage/KnownBlacklist.germline.bed
    touch output/sage/KnownBlacklist.germline.bed
    touch output/sage/KnownHotspots.somatic.vcf.gz
    touch output/sage/ActionableCodingPanel.bed.gz
    touch output/sage/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz
    touch output/sage/KnownBlacklist.germline.vcf.gz
    touch output/sage/KnownHotspots.somatic.vcf.gz
    touch output/sage/SageGermlinePon.98x.tsv.gz
    touch output/sage/clinvar.vcf.gz

    mkdir -p output/lilac/

    mkdir -p output/virusbreakend/

    mkdir -p output/virusinterpreter/
    touch output/virusinterpreter/taxonomy_db.tsv
    touch output/virusinterpreter/virus_reporting_db.tsv

    mkdir -p output/purple/
    touch output/purple/cohort_germline_del_freq.csv

    mkdir -p output/gene_panel/
    touch output/gene_panel/DriverGenePanel.tsv

    mkdir -p output/ensembl_data_cache/

    mkdir -p output/known_fusions/
    touch output/known_fusions/known_fusion_data.csv
    touch output/known_fusions/known_fusions.bedpe

    mkdir -p output/mappability/
    touch output/mappability/mappability_150.bed.gz
    """
}
