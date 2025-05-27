//
// Prepare reference data as required
//

import Constants

include { BWAMEM2_INDEX         } from '../../../modules/nf-core/bwamem2/index/main'
include { BWA_INDEX             } from '../../../modules/nf-core/bwa/index/main'
include { SAMTOOLS_DICT         } from '../../../modules/nf-core/samtools/dict/main'
include { SAMTOOLS_FAIDX        } from '../../../modules/nf-core/samtools/faidx/main'
include { STAR_GENOMEGENERATE   } from '../../../modules/nf-core/star/genomegenerate/main'

include { CUSTOM_EXTRACTTARBALL as DECOMP_BWAMEM2_INDEX    } from '../../../modules/local/custom/extract_tarball/main'
include { CUSTOM_EXTRACTTARBALL as DECOMP_GRIDSS_INDEX     } from '../../../modules/local/custom/extract_tarball/main'
include { CUSTOM_EXTRACTTARBALL as DECOMP_HMF_DATA         } from '../../../modules/local/custom/extract_tarball/main'
include { CUSTOM_EXTRACTTARBALL as DECOMP_PANEL_DATA       } from '../../../modules/local/custom/extract_tarball/main'
include { CUSTOM_EXTRACTTARBALL as DECOMP_STAR_INDEX       } from '../../../modules/local/custom/extract_tarball/main'
include { GATK4_BWA_INDEX_IMAGE                            } from '../../../modules/local/gatk4/bwaindeximage/main'
include { GRIDSS_INDEX                                     } from '../../../modules/local/gridss/index/main'
include { WRITE_REFERENCE_DATA                             } from '../../../modules/local/custom/write_reference_data/main'

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
    ch_genome_fasta = Channel.fromPath(params.ref_data_genome_fasta)
    ch_genome_version = Channel.value(params.genome_version)
    run_virusinterpreter = run_config.mode !== Constants.RunMode.TARGETED && run_config.stages.virusinterpreter

    //
    // Set .fai and .dict indexes, create if required
    //
    ch_genome_fai = getRefFileChannel('ref_data_genome_fai')
    if (!params.ref_data_genome_fai) {
        SAMTOOLS_FAIDX(ch_genome_fasta)
        ch_genome_fai = SAMTOOLS_FAIDX.out.fai
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    ch_genome_dict = getRefFileChannel('ref_data_genome_dict')
    if (!params.ref_data_genome_dict) {
        SAMTOOLS_DICT(ch_genome_fasta)
        ch_genome_dict = SAMTOOLS_DICT.out.dict
        ch_versions = ch_versions.mix(SAMTOOLS_DICT.out.versions)
    }

    ch_genome_img = getRefFileChannel('ref_data_genome_img')
    if (!params.ref_data_genome_img) {
        GATK4_BWA_INDEX_IMAGE(ch_genome_fasta)
        ch_genome_img = GATK4_BWA_INDEX_IMAGE.out.img
        ch_versions = ch_versions.mix(GATK4_BWA_INDEX_IMAGE.out.versions)
    }

    //
    // Set bwa-mem2 index, unpack or create if required
    //
    ch_genome_bwamem2_index = Channel.empty()
    if (run_config.has_dna_fastq && run_config.stages.alignment) {
        if (!params.ref_data_genome_bwamem2_index) {

            BWAMEM2_INDEX(
                ch_genome_fasta,
                params.ref_data_genome_alt ? file(params.ref_data_genome_alt) : [],
            )
            ch_genome_bwamem2_index = BWAMEM2_INDEX.out.index
            ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

        } else if (params.ref_data_genome_bwamem2_index.endsWith('.tar.gz')) {

            ch_genome_bwamem2_index_inputs = Channel.fromPath(params.ref_data_genome_bwamem2_index)
                .map { [[id: "${it.name.replaceAll('\\.tar\\.gz$', '')}"], it] }

            DECOMP_BWAMEM2_INDEX(ch_genome_bwamem2_index_inputs)
            ch_genome_bwamem2_index = DECOMP_BWAMEM2_INDEX.out.extracted_dir

        } else {

            ch_genome_bwamem2_index = getRefFileChannel('ref_data_genome_bwamem2_index')

        }
    }

    //
    // Set GRIDSS index, unpack or create if required
    //
    ch_genome_gridss_index = Channel.empty()
    if (run_config.has_dna && (run_config.stages.gridss || run_virusinterpreter)) {
        if (!params.ref_data_genome_gridss_index) {

            BWA_INDEX(
                ch_genome_fasta,
                params.ref_data_genome_alt ? file(params.ref_data_genome_alt) : [],
            )
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

            GRIDSS_INDEX(
                ch_genome_fasta,
                ch_genome_fai,
                ch_genome_dict,
                BWA_INDEX.out.index,
            )
            ch_genome_gridss_index = GRIDSS_INDEX.out.index
            ch_versions = ch_versions.mix(GRIDSS_INDEX.out.versions)

        } else if (params.ref_data_genome_gridss_index.endsWith('.tar.gz')) {

            ch_genome_gridss_index_inputs = Channel.fromPath(params.ref_data_genome_gridss_index)
                .map { [[id: "${it.name.replaceAll('\\.tar\\.gz$', '')}"], it] }

            DECOMP_GRIDSS_INDEX(ch_genome_gridss_index_inputs)
            ch_genome_gridss_index = DECOMP_GRIDSS_INDEX.out.extracted_dir

        } else {

            ch_genome_gridss_index = getRefFileChannel('ref_data_genome_gridss_index')

        }
    }

    //
    // Set STAR index , unpack or create if required
    //
    ch_genome_star_index = Channel.empty()
    if (run_config.has_rna_fastq && run_config.stages.alignment) {
        if (!params.ref_data_genome_star_index) {

            STAR_GENOMEGENERATE(
                ch_genome_fasta,
                file(params.ref_data_genome_gtf),
            )
            ch_genome_star_index = STAR_GENOMEGENERATE.out.index
            ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

        } else if (params.ref_data_genome_star_index.endsWith('.tar.gz')) {

            ch_genome_star_index_inputs = Channel.fromPath(params.ref_data_genome_star_index)
                .map { [[id: "${it.name.replaceAll('\\.tar\\.gz$', '')}"], it] }

            DECOMP_STAR_INDEX(ch_genome_star_index_inputs)
            ch_genome_star_index = DECOMP_STAR_INDEX.out.extracted_dir

        } else {

            ch_genome_star_index = getRefFileChannel('ref_data_genome_star_index')

        }
    }

    //
    // Set HMF reference data, unpack if required
    //
    ch_hmf_data = Channel.empty()
    hmf_data_paths = params.hmf_data_paths[params.genome_version.toString()]
    if (params.ref_data_hmf_data_path.endsWith('tar.gz')) {

        ch_hmf_data_inputs = Channel.fromPath(params.ref_data_hmf_data_path)
            .map { [[id: "${it.name.replaceAll('\\.tar\\.gz$', '')}"], it] }

        DECOMP_HMF_DATA(ch_hmf_data_inputs)

        ch_hmf_data = DECOMP_HMF_DATA.out.extracted_dir
            .collect()
            .map { dir_list ->
                assert dir_list.size() == 1
                def dirpath = dir_list[0].toUriString()
                return createDataMap(hmf_data_paths, dirpath)
            }

    } else {

        ch_hmf_data = Channel.value(createDataMap(hmf_data_paths, params.ref_data_hmf_data_path))

    }

    //
    // Set panel reference data, unpack if required
    //
    ch_panel_data = Channel.empty()
    if (run_config.mode === Constants.RunMode.TARGETED) {

        panel_data_paths_versions = params.panel_data_paths[params.panel]
        panel_data_paths = panel_data_paths_versions[params.genome_version.toString()]

        if (params.ref_data_panel_data_path.endsWith('tar.gz')) {

            ch_panel_data_inputs = Channel.fromPath(params.ref_data_panel_data_path)
                .map { [[id: "${it.name.replaceAll('\\.tar\\.gz$', '')}"], it] }

            DECOMP_PANEL_DATA(ch_panel_data_inputs)

            ch_panel_data = DECOMP_PANEL_DATA.out.extracted_dir
                .collect()
                .map { dir_list ->
                    assert dir_list.size() == 1
                    def dirpath = dir_list[0].toUriString()
                    return createDataMap(panel_data_paths, dirpath)
                }

        } else {

            ch_panel_data = Channel.value(createDataMap(panel_data_paths, params.ref_data_panel_data_path))

        }
    }

    //
    // Write prepared reference data if requested
    //
    if (params.prepare_reference_only) {

        // Create channel of data files to stage (if not already local) and write
        ch_refdata = Channel.empty()
            .mix(
                ch_genome_fasta,
                ch_genome_fai,
                ch_genome_dict,
                ch_genome_img,
                ch_genome_bwamem2_index,
                ch_genome_gridss_index,
                ch_genome_star_index,
                // Also include base paths for hmf_data and panel_data
                Channel.empty()
                    .mix(
                        ch_hmf_data,
                        ch_panel_data,
                    )
                    .map { getDataBaseDirectory(it) }
            )

        WRITE_REFERENCE_DATA(
            ch_refdata,
            workflow.manifest.version,
        )

        // Clear all stages to prevent running any analysis
        run_config.stages = [:]
    }

    emit:
    genome_fasta         = ch_genome_fasta.first()         // path: genome_fasta
    genome_fai           = ch_genome_fai.first()           // path: genome_fai
    genome_dict          = ch_genome_dict.first()          // path: genome_dict
    genome_img           = ch_genome_img.first()           // path: genome_img
    genome_bwamem2_index = ch_genome_bwamem2_index.first() // path: genome_bwa-mem2_index
    genome_gridss_index  = ch_genome_gridss_index.first()  // path: genome_gridss_index
    genome_star_index    = ch_genome_star_index.first()    // path: genome_star_index
    genome_version       = ch_genome_version               // val:  genome_version

    hmf_data             = ch_hmf_data                     // map:  HMF data paths
    panel_data           = ch_panel_data                   // map:  Panel data paths

    versions             = ch_versions                     // channel: [ versions.yml ]
}

def getRefFileChannel(key) {
    def fp = params.get(key) ? file(params.getAt(key)) : []
    return Channel.of(fp)
}

def createDataMap(entries, ref_data_path) {
    return entries
        .collectEntries { name, path ->
            def ref_data_file = path == [] ? [] : getRefdataFile(path, ref_data_path)
            return [name, ref_data_file]
        }
}

def getRefdataFile(filepath, ref_data_path) {
    def data_path_noslash = ref_data_path.toString().replaceAll('/$', '')
    return file("${data_path_noslash}/${filepath}", checkIfExists: true)
}

def getDataBaseDirectory(data) {
    def c = []
    data
        .findAll { it.value }
        .collect { it.value.toUriString().getChars() }
        .transpose()
        .findIndexOf {
            def cs = it.unique()
            if (cs.size() != 1) return true
            c << cs.pop()
            return false
        }
    return file("${c.join('')}")
}
