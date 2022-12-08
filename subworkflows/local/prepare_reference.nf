//
// Prepare reference data as required
//

include { SAMTOOLS_FAIDX         } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_DICT          } from '../../modules/nf-core/samtools/dict/main'
include { BWA_INDEX              } from '../../modules/nf-core/bwa/index/main'

include { EXTRACT_TARBALL as EXTRACT_TARBALL_HMF_BUNDLE      } from '../../modules/local/custom/extract_tarball/main'
include { EXTRACT_TARBALL as EXTRACT_TARBALL_VIRUSBREAKENDDB } from '../../modules/local/custom/extract_tarball/main'
include { INDEX as GRIDSS_BWA_INDEX_IMAGE                    } from '../../modules/local/gridss/index/main'
include { INDEX as GRIDSS_INDEX                              } from '../../modules/local/gridss/index/main'

workflow PREPARE_REFERENCE {
    take:
        run // map: stages to run

    main:
        // Channel for version.yml files
        ch_versions = Channel.empty()

        //
        // Set reference genome FASTA for consistency
        //
        ch_genome_fasta = file(params.ref_data_genome_fasta)

        //
        // Create .fai and .dict for reference genome if required
        //
        // The fai and dict files should always be present if using a genome preset. These are
        // always created where they are not present without checking processes to run given they
        // are used in numerous processes and have a neglibile cost to generate.
        ch_genome_fai = file(params.ref_data_genome_fai)
        ch_genome_dict = file(params.ref_data_genome_dict)
        if (!params.ref_data_genome_fai) {
            SAMTOOLS_FAIDX([:], ch_genome_fasta)
            ch_genome_fai = SAMTOOLS_FAIDX.out.fai
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        }
        if (!params.ref_data_genome_dict) {
            SAMTOOLS_DICT([:], ch_genome_fasta)
            ch_genome_dict = SAMTOOLS_DICT.out.dict
            ch_versions = ch_versions.mix(SAMTOOLS_DICT.out.versions)
        }

        //
        // Create BWA index, BWA index image, and GRIDSS index for reference genome if required
        //
        ch_genome_bwa_index = file(params.ref_data_genome_bwa_index)
        ch_genome_bwa_index_image = file(params.ref_data_genome_bwa_index_image)
        ch_genome_gridss_index = file(params.ref_data_genome_gridss_index)
        if (run.gridss || run.virusinterpreter) {
            if (!params.ref_data_genome_bwa_index) {
                BWA_INDEX([[:], ch_genome_fasta])
                ch_genome_bwa_index = BWA_INDEX.out.index.map { it[1] }
                ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
            }
            if (!params.ref_data_genome_bwa_index_image) {
                GRIDSS_BWA_INDEX_IMAGE(
                    ch_genome_fasta,
                    ch_genome_fai,
                    ch_genome_dict,
                    ch_genome_bwa_index,
                    [],
                    ['bwa_index_image'],
                )
                ch_genome_bwa_index_image = GRIDSS_BWA_INDEX_IMAGE.out.img
                ch_versions = ch_versions.mix(GRIDSS_BWA_INDEX_IMAGE.out.versions)
            }
            if (!params.ref_data_genome_gridss_index) {
                GRIDSS_INDEX(
                    ch_genome_fasta,
                    ch_genome_fai,
                    ch_genome_dict,
                    ch_genome_bwa_index,
                    ch_genome_bwa_index_image,
                    ['gridss_index'],
                )
                ch_genome_gridss_index = GRIDSS_INDEX.out.index
                ch_versions = ch_versions.mix(GRIDSS_INDEX.out.versions)
            }
        }

        //
        // Set VIRUSBreakend database paths / stage, unpack if required
        //
        ch_virusbreakenddb = Channel.empty()
        if (run.virusinterpreter) {
            if (params.ref_data_virusbreakenddb_path.endsWith('.tar.gz')) {
                ch_virusbreakenddb_inputs = [
                    [id: 'virusbreakenddb'],
                    file(params.ref_data_virusbreakenddb_path),
                ]
                EXTRACT_TARBALL_VIRUSBREAKENDDB(ch_virusbreakenddb_inputs)
                ch_virusbreakenddb = EXTRACT_TARBALL_VIRUSBREAKENDDB.out.dir
            } else {
                ch_virusbreakenddb = file(params.ref_data_virusbreakenddb_path)
            }
        }

        //
        // Set HMF reference paths / stage, unpack if required
        //
        // NOTE(SW): requiring all HMF reference data for now
        if (params.ref_data_hmf_bundle.endsWith('.tar.gz')) {
            // Decompress and set paths
            EXTRACT_TARBALL_HMF_BUNDLE([[id: 'hmf_bundle'], file(params.ref_data_hmf_bundle)])

            // Obtain paths and convert queue channel to value channel
            // NOTE(SW): any relevant HMF reference file parameter explicitly set in config will take priority
            ch_hmf_data = EXTRACT_TARBALL_HMF_BUNDLE.out.dir
                .collect()
                .map { dir_list ->
                    assert dir_list.size() == 1
                    return createHmfDataMap(dir_list[0], false /* params_only */)
                }
        } else if (params.ref_data_hmf_bundle) {
            // If provided as path to directory, set paths
            ch_hmf_data = createHmfDataMap(params.ref_data_hmf_bundle, false /* params_only */)
        } else {
            // If no HMF data bundle is supplied we construct from *only* params
            ch_hmf_data = createHmfDataMap(null, true /* params_only */)
        }

    emit:
        genome_fasta           = ch_genome_fasta                // path: genome_fasta
        genome_fai             = ch_genome_fai                  // path: genome_fai
        genome_dict            = ch_genome_dict                 // path: genome_dict
        genome_bwa_index       = ch_genome_bwa_index            // path: genome_bwa_index
        genome_bwa_index_image = ch_genome_bwa_index_image      // path: genome_bwa_index_image
        genome_gridss_index    = ch_genome_gridss_index         // path: genome_gridss_index
        genome_version         = params.ref_data_genome_version // val:  genome_version

        virusbreakenddb        = ch_virusbreakenddb             // path: VIRUSBreakend database
        hmf_data               = ch_hmf_data                    // map:  HMF data paths

        versions               = ch_versions                    // channel: [versions.yml]
}

def createHmfDataMap(hmf_bundle_dir, params_only) {
    return [
        // AMBER
        'amber_loci':                   getHmfDataFileObject('ref_data_amber_loci',                   'AMBER_LOCI',                   hmf_bundle_dir, params_only),
        // COBALT
        'cobalt_gc_profile':            getHmfDataFileObject('ref_data_cobalt_gc_profile',            'COBALT_GC_PROFILE',            hmf_bundle_dir, params_only),
        // CUPPA
        'cuppa':                        getHmfDataFileObject('ref_data_cuppa',                        'CUPPA',                        hmf_bundle_dir, params_only),
        // SVPREP
        'sv_prep_blacklist':            getHmfDataFileObject('ref_data_sv_prep_blacklist',            'SV_PREP_BLACKLIST',            hmf_bundle_dir, params_only),
        // GRIDSS, GRIPSS
        'gridss_blacklist':             getHmfDataFileObject('ref_data_gridss_blacklist',             'GRIDSS_BLACKLIST',             hmf_bundle_dir, params_only),
        'gridss_breakend_pon':          getHmfDataFileObject('ref_data_gridss_breakend_pon',          'GRIDSS_BREAKEND_PON',          hmf_bundle_dir, params_only),
        'gridss_breakpoint_pon':        getHmfDataFileObject('ref_data_gridss_breakpoint_pon',        'GRIDSS_BREAKPOINT_PON',        hmf_bundle_dir, params_only),
        'repeat_masker_file':           getHmfDataFileObject('ref_data_repeat_masker_file',           'REPEAT_MASKER_FILE',           hmf_bundle_dir, params_only),
        // Isofox
        'isofox_exp_counts':            getHmfDataFileObject('ref_data_isofox_exp_counts',            'ISOFOX_EXP_COUNTS',            hmf_bundle_dir, params_only),
        'isofox_exp_gc_ratios':         getHmfDataFileObject('ref_data_isofox_exp_gc_ratios',         'ISOFOX_EXP_GC_RATIOS',         hmf_bundle_dir, params_only),
        // LINX
        'linx_fragile_sites':           getHmfDataFileObject('ref_data_linx_fragile_sites',           'LINX_FRAGILE_SITES',           hmf_bundle_dir, params_only),
        'linx_lines':                   getHmfDataFileObject('ref_data_linx_lines',                   'LINX_LINES',                   hmf_bundle_dir, params_only),
        // SAGE, PAVE
        'sage_blacklist_bed':           getHmfDataFileObject('ref_data_sage_blacklist_bed',           'SAGE_BLACKLIST_BED',           hmf_bundle_dir, params_only),
        'sage_blacklist_vcf':           getHmfDataFileObject('ref_data_sage_blacklist_vcf',           'SAGE_BLACKLIST_VCF',           hmf_bundle_dir, params_only),
        'sage_coding_panel':            getHmfDataFileObject('ref_data_sage_coding_panel',            'SAGE_CODING_PANEL',            hmf_bundle_dir, params_only),
        'sage_high_confidence':         getHmfDataFileObject('ref_data_sage_high_confidence',         'SAGE_HIGH_CONFIDENCE',         hmf_bundle_dir, params_only),
        'sage_known_hotspots_germline': getHmfDataFileObject('ref_data_sage_known_hotspots_germline', 'SAGE_KNOWN_HOTSPOTS_GERMLINE', hmf_bundle_dir, params_only),
        'sage_known_hotspots_somatic':  getHmfDataFileObject('ref_data_sage_known_hotspots_somatic',  'SAGE_KNOWN_HOTSPOTS_SOMATIC',  hmf_bundle_dir, params_only),
        'sage_pon_file':                getHmfDataFileObject('ref_data_sage_pon_file',                'SAGE_PON_FILE',                hmf_bundle_dir, params_only),
        'clinvar_vcf':                  getHmfDataFileObject('ref_data_clinvar_vcf',                  'CLINVAR_VCF',                  hmf_bundle_dir, params_only),
        // SIGS
        'sigs_signatures':              getHmfDataFileObject('ref_data_sigs_signatures',              'SIGS_SIGNATURES',              hmf_bundle_dir, params_only),
        // LILAC
        'lilac_resource_dir':           getHmfDataFileObject('ref_data_lilac_resource_dir',           'LILAC_RESOURCE_DIR',           hmf_bundle_dir, params_only),
        // Virus Interpreter
        'virus_taxonomy':               getHmfDataFileObject('ref_data_virus_taxonomy',               'VIRUSINTERPRETER_TAXONOMY',    hmf_bundle_dir, params_only),
        'virus_reporting':              getHmfDataFileObject('ref_data_virus_reporting',              'VIRUSINTERPRETER_REPORTING',   hmf_bundle_dir, params_only),
        // Misc
        'purple_germline_del':          getHmfDataFileObject('ref_data_purple_germline_del',          'PURPLE_GERMLINE_DEL',          hmf_bundle_dir, params_only),
        'driver_gene_panel':            getHmfDataFileObject('ref_data_driver_gene_panel',            'DRIVER_GENE_PANEL',            hmf_bundle_dir, params_only),
        'ensembl_data_dir':             getHmfDataFileObject('ref_data_ensembl_data_dir',             'ENSEMBL_DATA_DIR',             hmf_bundle_dir, params_only),
        'known_fusion_data':            getHmfDataFileObject('ref_data_known_fusion_data',            'KNOWN_FUSION_DATA',            hmf_bundle_dir, params_only),
        'known_fusions':                getHmfDataFileObject('ref_data_known_fusions',                'KNOWN_FUSIONS',                hmf_bundle_dir, params_only),
        'mappability_bed':              getHmfDataFileObject('ref_data_mappability_bed',              'MAPPABILITY_BED',              hmf_bundle_dir, params_only),
    ]
}

def getHmfDataFileObject(pk, hk, base_dir, params_only) {
    if (params.containsKey(pk)) {
        return file(params.getAt(pk), checkIfExists: true)
    } else if (params_only) {
        assert false : "TODO(SW): more helpful message about missing param ${pk}"
    } else if (!Constants.HMF_DATA_PATHS.containsKey(hk)) {
        assert false : "bad key for Constants.HMF_DATA_PATHS ${hk}"
    } else {
      def base_dir_noslash = base_dir.toString().replaceAll('/$', '')
      return file("${base_dir_noslash}/${Constants.HMF_DATA_PATHS.getAt(hk)}", checkIfExists: true)
    }
}
