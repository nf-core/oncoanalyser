//
// Prepare reference data as required
//

// NOTE(SW): BWA MEM2 indexes are required and are not created
// TODO(SW): consider removing index creation since it's unlikely to be used, replace with documentation

import Constants

include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_DICT  } from '../../../modules/nf-core/samtools/dict/main'
include { BWA_INDEX      } from '../../../modules/nf-core/bwa/index/main'

include { CUSTOM_EXTRACTTARBALL as DECOMP_BWA_INDEX        } from '../../../modules/local/custom/extract_tarball/main'
include { CUSTOM_EXTRACTTARBALL as DECOMP_HMF_DATA         } from '../../../modules/local/custom/extract_tarball/main'
include { CUSTOM_EXTRACTTARBALL as DECOMP_PANEL_DATA       } from '../../../modules/local/custom/extract_tarball/main'
include { CUSTOM_EXTRACTTARBALL as DECOMP_STAR_INDEX       } from '../../../modules/local/custom/extract_tarball/main'
include { CUSTOM_EXTRACTTARBALL as DECOMP_VIRUSBREAKEND_DB } from '../../../modules/local/custom/extract_tarball/main'
include { GRIDSS_INDEX as GRIDSS_BWA_INDEX_IMAGE           } from '../../../modules/local/gridss/index/main'
include { GRIDSS_INDEX as GRIDSS_INDEX                     } from '../../../modules/local/gridss/index/main'

workflow PREPARE_REFERENCE {
    take:
        run_config // channel: [mandatory] run configuration

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        //
        // Set some variables for brevity
        //
        ch_genome_fasta = file(params.ref_data.genome_fasta)
        run_virusinterpreter = run_config.mode !== Constants.RunMode.TARGETED && run_config.stages.virusinterpreter

        //
        // Create .fai and .dict for reference genome if required
        //
        // The fai and dict files should always be present if using a genome preset. These are
        // always created where they are not present without checking processes to run given they
        // are used in numerous processes and have a neglibile cost to generate.
        ch_genome_fai = params.ref_data.genome_fai ? file(params.ref_data.genome_fai) : []
        ch_genome_dict = params.ref_data.genome_dict ? file(params.ref_data.genome_dict) : []
        if (!params.ref_data.genome_fai) {
            SAMTOOLS_FAIDX([[:], ch_genome_fasta])
            ch_genome_fai = SAMTOOLS_FAIDX.out.fai.map { meta, fai -> fai }
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        }
        if (!params.ref_data.genome_dict) {
            SAMTOOLS_DICT([[:], ch_genome_fasta])
            ch_genome_dict = SAMTOOLS_DICT.out.dict.map { meta, dict -> dict }
            ch_versions = ch_versions.mix(SAMTOOLS_DICT.out.versions)
        }

        //
        // Create BWA index, BWA index image, and GRIDSS index for reference genome if required
        //
        ch_genome_bwa_index = params.ref_data.genome_bwa_index ? file(params.ref_data.genome_bwa_index) : []
        ch_genome_bwa_index_image = params.ref_data.genome_gridss_index ? file(params.ref_data.genome_bwa_index_image) : []
        ch_genome_gridss_index = params.ref_data.genome_gridss_index ? file(params.ref_data.genome_gridss_index) : []
        if (run_config.has_dna && (run_config.stages.gridss || run_virusinterpreter)) {
            // NOTE(SW): the BWA index directory can be provided as a compressed tarball
            if (!params.ref_data.genome_bwa_index) {
                BWA_INDEX([[:], ch_genome_fasta])
                ch_genome_bwa_index = BWA_INDEX.out.index.map { meta, index -> index }
                ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
            } else if (params.ref_data.genome_bwa_index.endsWith('.tar.gz')) {
                ch_genome_bwa_index_inputs = [
                    [id: 'bwa_index'],
                    file(params.ref_data.genome_bwa_index),
                ]
                DECOMP_BWA_INDEX(ch_genome_bwa_index_inputs)
                ch_genome_bwa_index = DECOMP_BWA_INDEX.out.dir
            }

            if (!params.ref_data.genome_bwa_index_image) {
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
            if (!params.ref_data.genome_gridss_index) {
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

        // Explicitly set BWA MEM2 index file inputs
        ch_genome_bwa_index_bseq = file(params.ref_data.genome_bwa_index_bseq)
        ch_genome_bwa_index_biidx = file(params.ref_data.genome_bwa_index_biidx)

        //
        // Decompress STAR index
        //
        ch_genome_star_index = params.ref_data.genome_star_index ? file(params.ref_data.genome_star_index) : []
        if (run_config.has_rna_fastq && run_config.stages.alignment && params.ref_data.genome_star_index.endsWith('.tar.gz')) {
                ch_genome_star_index_inputs = [
                    [id: 'star_index'],
                    file(params.ref_data.genome_star_index),
                ]
                DECOMP_STAR_INDEX(ch_genome_star_index_inputs)
                ch_genome_star_index = DECOMP_STAR_INDEX.out.dir
        }

        //
        // Set VIRUSBreakend database path / stage, unpack if required
        //
        ch_virusbreakenddb = Channel.empty()
        if (run_config.has_dna && run_virusinterpreter) {
            if (params.ref_data.virusbreakenddb_path.endsWith('.tar.gz')) {
                ch_virusbreakenddb_inputs = [
                    [id: 'virusbreakenddb'],
                    file(params.ref_data.virusbreakenddb_path),
                ]
                DECOMP_VIRUSBREAKEND_DB(ch_virusbreakenddb_inputs)
                ch_virusbreakenddb = DECOMP_VIRUSBREAKEND_DB.out.dir
            } else {
                ch_virusbreakenddb = file(params.ref_data.virusbreakenddb_path)
            }
        }

        //
        // Set HMF reference paths / stage, unpack if required
        //
        hmf_data_paths = params.hmf_data_paths[params.ref_data.genome_version]
        if (params.ref_data.hmf_data_path.endsWith('tar.gz')) {
            ch_hmf_data_inputs = [
                [id: 'hmf_data'],
                file(params.ref_data.hmf_data_path),
            ]
            DECOMP_HMF_DATA(ch_hmf_data_inputs)

            ch_hmf_data = DECOMP_HMF_DATA.out.dir
                .collect()
                .map { dir_list ->
                    assert dir_list.size() == 1
                    def dirpath = dir_list[0].toUriString()
                    return createDataMap(hmf_data_paths, dirpath)
                }
        } else {
            ch_hmf_data = createDataMap(hmf_data_paths, params.ref_data.hmf_data_path)
        }

        //
        // Set panel reference paths / stage, unpack if required
        //
        ch_panel_data = Channel.empty()
        if (run_config.mode === Constants.RunMode.TARGETED) {

            // NOTE(SW): consider approach to implement custom panel support

            panel_data_paths_versions = params.panel_data_paths[params.panel]
            panel_data_paths = panel_data_paths_versions[params.ref_data.genome_version]

            if (params.ref_data.panel_data_path.endsWith('tar.gz')) {
                ch_panel_data_inputs = [
                    [id: 'panel_data'],
                    file(params.ref_data.panel_data_path),
                ]
                DECOMP_PANEL_DATA(ch_panel_data_inputs)

                ch_panel_data = DECOMP_PANEL_DATA.out.dir
                    .collect()
                    .map { dir_list ->
                        assert dir_list.size() == 1
                        def dirpath = dir_list[0].toUriString()
                        return createDataMap(panel_data_paths, dirpath)
                    }
            } else {
                ch_panel_data = createDataMap(panel_data_paths, params.ref_data.panel_data_path)
            }
        }

    emit:
        genome_fasta           = ch_genome_fasta                // path: genome_fasta
        genome_fai             = ch_genome_fai                  // path: genome_fai
        genome_dict            = ch_genome_dict                 // path: genome_dict
        genome_bwa_index       = ch_genome_bwa_index            // path: genome_bwa_index
        genome_bwa_index_bseq  = ch_genome_bwa_index_bseq       // path: genome_bwa_index_bseq
        genome_bwa_index_biidx = ch_genome_bwa_index_biidx      // path: genome_bwa_index_biidx
        genome_bwa_index_image = ch_genome_bwa_index_image      // path: genome_bwa_index_image
        genome_gridss_index    = ch_genome_gridss_index         // path: genome_gridss_index
        genome_star_index      = ch_genome_star_index           // path: genome_star_index
        genome_version         = params.ref_data.genome_version // val:  genome_version

        virusbreakenddb        = ch_virusbreakenddb             // path: VIRUSBreakend database
        hmf_data               = ch_hmf_data                    // map:  HMF data paths
        panel_data             = ch_panel_data                  // map:  Panel data paths

        versions               = ch_versions                    // channel: [ versions.yml ]
}

def createDataMap(entries, ref_data_path) {
    return entries
        .collectEntries { name, path ->
            def ref_data_file = getRefdataFile(path, ref_data_path)
            return [name, ref_data_file]
        }
}

def getRefdataFile(filepath, ref_data_path) {
    def data_path_noslash = ref_data_path.toString().replaceAll('/$', '')
    return file("${data_path_noslash}/${filepath}", checkIfExists: true)
}
