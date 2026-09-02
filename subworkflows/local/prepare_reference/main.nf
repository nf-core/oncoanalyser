//
// Prepare reference data as required
//

include { BWAMEM2_INDEX         } from '../../../modules/nf-core/bwamem2/index/main'
include { BWA_INDEX             } from '../../../modules/nf-core/bwa/index/main'
include { SAMTOOLS_DICT         } from '../../../modules/nf-core/samtools/dict/main'
include { SAMTOOLS_FAIDX        } from '../../../modules/nf-core/samtools/faidx/main'
include { GATK4_BWA_INDEX_IMAGE } from '../../../modules/local/gatk4/bwaindeximage/main'
include { STAR_GENOMEGENERATE   } from '../../../modules/nf-core/star/genomegenerate/main'
include { GRIDSS_INDEX          } from '../../../modules/local/gridss/index/main'

include { CUSTOM_EXTRACTTARBALL as DECOMP_BWAMEM2_INDEX } from '../../../modules/local/custom/extract_tarball/main'
include { CUSTOM_EXTRACTTARBALL as DECOMP_GRIDSS_INDEX  } from '../../../modules/local/custom/extract_tarball/main'
include { CUSTOM_EXTRACTTARBALL as DECOMP_HMF_DATA      } from '../../../modules/local/custom/extract_tarball/main'
include { CUSTOM_EXTRACTTARBALL as DECOMP_PANEL_DATA    } from '../../../modules/local/custom/extract_tarball/main'
include { CUSTOM_EXTRACTTARBALL as DECOMP_STAR_INDEX    } from '../../../modules/local/custom/extract_tarball/main'

include { WRITE_REFERENCE_DATA as WRITE_FASTA           } from '../../../modules/local/custom/write_reference_data/main'
include { WRITE_REFERENCE_DATA as WRITE_FAI             } from '../../../modules/local/custom/write_reference_data/main'
include { WRITE_REFERENCE_DATA as WRITE_DICT            } from '../../../modules/local/custom/write_reference_data/main'
include { WRITE_REFERENCE_DATA as WRITE_IMG             } from '../../../modules/local/custom/write_reference_data/main'
include { WRITE_REFERENCE_DATA as WRITE_BWA_INDEX       } from '../../../modules/local/custom/write_reference_data/main'
include { WRITE_REFERENCE_DATA as WRITE_GRIDSS_INDEX    } from '../../../modules/local/custom/write_reference_data/main'
include { WRITE_REFERENCE_DATA as WRITE_STAR_INDEX      } from '../../../modules/local/custom/write_reference_data/main'
include { WRITE_REFERENCE_DATA as WRITE_HMF_DATA        } from '../../../modules/local/custom/write_reference_data/main'
include { WRITE_REFERENCE_DATA as WRITE_PANEL_DATA      } from '../../../modules/local/custom/write_reference_data/main'

workflow PREPARE_REFERENCE {
    take:
    prep_config // channel: [mandatory] configuration indicating which reference data is required
    run_config
    params

    main:
    //
    // Set .fasta and main genome indexes, create if required
    //
    ch_genome_version = channel.value(params.genome_version)

    ch_genome_fasta = channel.empty()
    if (prep_config.require_fasta) {
        ch_genome_fasta = channel.fromPath(params.ref_data_genome_fasta)
    }

    ch_genome_fai = channel.empty()
    if (prep_config.require_fai) {

        if (! params.ref_data_genome_fai) {
            SAMTOOLS_FAIDX(ch_genome_fasta)
            ch_genome_fai = channel.topic('samtools_fai')
        } else {
            ch_genome_fai = channel.fromPath(params.ref_data_genome_fai)
        }
    }

    ch_genome_dict = channel.empty()
    if (prep_config.require_dict) {

        if (! params.ref_data_genome_dict) {
            SAMTOOLS_DICT(ch_genome_fasta)
            ch_genome_dict = channel.topic('samtools_dict')
        } else {
            ch_genome_dict = channel.fromPath(params.ref_data_genome_dict)
        }
    }

    ch_genome_img = channel.empty()
    if (prep_config.require_img) {

        if (! params.ref_data_genome_img) {
            GATK4_BWA_INDEX_IMAGE(ch_genome_fasta)
            ch_genome_img = channel.topic('gatk4_bwa_index_img')
        } else {
            ch_genome_img = channel.fromPath(params.ref_data_genome_img)
        }
    }

    //
    // Set bwa-mem2 index, unpack or create if required
    //
    ch_genome_bwamem2_index = channel.empty()
    if (prep_config.require_bwamem2_index) {

        if (! params.ref_data_genome_bwamem2_index) {

            BWAMEM2_INDEX(
                ch_genome_fasta,
                params.ref_data_genome_alt ? file(params.ref_data_genome_alt) : [],
            )
            ch_genome_bwamem2_index = channel.topic('bwamem2_index')

        } else if (params.ref_data_genome_bwamem2_index.endsWith('.tar.gz')) {

            ch_genome_bwamem2_index_inputs = channel.of(params.ref_data_genome_bwamem2_index)
                .map { fp_str -> def fp = file(fp_str); return [[topic_key: fp_str, id: "${fp.name.replaceAll('\\.tar\\.gz\$', '')}"], fp] }

            DECOMP_BWAMEM2_INDEX(ch_genome_bwamem2_index_inputs)
            ch_genome_bwamem2_index = channel.topic('extracted_dir')
                .filter { meta, _dir -> meta.topic_key == params.ref_data_genome_bwamem2_index }
                .map { _meta, dir -> dir }

        } else {

            ch_genome_bwamem2_index = channel.fromPath(params.ref_data_genome_bwamem2_index)

        }
    }

    //
    // Set GRIDSS index, unpack or create if required
    //
    ch_genome_gridss_index = channel.empty()
    if (prep_config.require_gridss_index) {

        if (! params.ref_data_genome_gridss_index) {

            BWA_INDEX(
                ch_genome_fasta,
                params.ref_data_genome_alt ? file(params.ref_data_genome_alt) : [],
            )

            GRIDSS_INDEX(
                ch_genome_fasta,
                ch_genome_fai,
                ch_genome_dict,
                channel.topic('bwa_index'),
            )
            ch_genome_gridss_index = channel.topic('gridss_index')

        } else if (params.ref_data_genome_gridss_index.endsWith('.tar.gz')) {

            ch_genome_gridss_index_inputs = channel.of(params.ref_data_genome_gridss_index)
                .map { fp_str -> def fp = file(fp_str); return [[topic_key: fp_str, id: "${fp.name.replaceAll('\\.tar\\.gz\$', '')}"], fp] }

            DECOMP_GRIDSS_INDEX(ch_genome_gridss_index_inputs)
            ch_genome_gridss_index = channel.topic('extracted_dir')
                .filter { meta, _dir -> meta.topic_key == params.ref_data_genome_gridss_index }
                .map { _meta, dir -> dir }

        } else {

            ch_genome_gridss_index = channel.fromPath(params.ref_data_genome_gridss_index)

        }
    }

    //
    // Set STAR index, unpack or create if required
    //
    ch_genome_star_index = channel.empty()
    if (prep_config.require_star_index) {

        if (! params.ref_data_genome_star_index) {

            STAR_GENOMEGENERATE(
                ch_genome_fasta,
                file(params.ref_data_genome_gtf),
            )
            ch_genome_star_index = channel.topic('star_index')

        } else if (params.ref_data_genome_star_index.endsWith('.tar.gz')) {

            ch_genome_star_index_inputs = channel.of(params.ref_data_genome_star_index)
                .map { fp_str -> def fp = file(fp_str); return [[topic_key: fp_str, id: "${fp.name.replaceAll('\\.tar\\.gz\$', '')}"], fp] }

            DECOMP_STAR_INDEX(ch_genome_star_index_inputs)
            ch_genome_star_index = channel.topic('extracted_dir')
                .filter { meta, _dir -> meta.topic_key == params.ref_data_genome_star_index }
                .map { _meta, dir -> dir }

        } else {

            ch_genome_star_index = channel.fromPath(params.ref_data_genome_star_index)

        }
    }

    //
    // Set HMF reference data, unpack if required
    //
    ch_hmf_data = channel.empty()
    if (prep_config.require_hmftools_data) {

        hmf_data_paths = params.hmf_data_paths[params.genome_version.toString()]

        if (params.ref_data_hmf_data_path.endsWith('tar.gz')) {

            ch_hmf_data_inputs = channel.of(params.ref_data_hmf_data_path)
                .map { fp_str -> def fp = file(fp_str); return [[topic_key: fp_str, id: "${fp.name.replaceAll('\\.tar\\.gz\$', '')}"], fp] }

            DECOMP_HMF_DATA(ch_hmf_data_inputs)

            ch_hmf_data = channel.topic('extracted_dir')
                .filter { meta, _dir -> meta.topic_key == params.ref_data_hmf_data_path }
                .map { _meta, dir -> dir }
                .collect()
                .map { dir_list ->
                    assert dir_list.size() == 1
                    def dirpath = dir_list[0].toUriString()
                    return createDataMap(hmf_data_paths, dirpath)
                }

        } else {

            ch_hmf_data = channel.value(createDataMap(hmf_data_paths, params.ref_data_hmf_data_path))

        }

    }

    //
    // Set panel reference data, unpack if required
    //
    ch_panel_data = channel.empty()
    if (prep_config.require_panel_data) {

        panel_data_paths_versions = params.panel_data_paths[params.panel.toLowerCase()]
        panel_data_paths = panel_data_paths_versions[params.genome_version.toString()]

        if (params.ref_data_panel_data_path.endsWith('tar.gz')) {

            ch_panel_data_inputs = channel.of(params.ref_data_panel_data_path)
                .map { fp_str -> def fp = file(fp_str); return [[topic_key: fp_str, id: "${fp.name.replaceAll('\\.tar\\.gz\$', '')}"], fp] }

            DECOMP_PANEL_DATA(ch_panel_data_inputs)

            ch_panel_data = channel.topic('extracted_dir')
                .filter { meta, dir -> meta.topic_key == params.ref_data_panel_data_path }
                .map { meta, dir -> dir }
                .collect()
                .map { dir_list ->
                    assert dir_list.size() == 1
                    def dirpath = dir_list[0].toUriString()
                    return createDataMap(panel_data_paths, dirpath)
                }

        } else {

            ch_panel_data = channel.value(createDataMap(panel_data_paths, params.ref_data_panel_data_path))

        }
    }

    //
    // Write prepared reference data if requested
    //
    if (prep_config.prepare_ref_data_only || params.prepare_reference_only) {

        WRITE_FASTA(ch_genome_fasta)
        WRITE_FAI(ch_genome_fai)
        WRITE_DICT(ch_genome_dict)
        WRITE_IMG(ch_genome_img)
        WRITE_BWA_INDEX(ch_genome_bwamem2_index)
        WRITE_GRIDSS_INDEX(ch_genome_gridss_index)
        WRITE_STAR_INDEX(ch_genome_star_index)

        WRITE_HMF_DATA(ch_hmf_data.map { f -> getDataBaseDirectory(f) })
        WRITE_PANEL_DATA(ch_panel_data.map { f -> getDataBaseDirectory(f) })

        // Clear all stages to prevent running any analysis when driving by samplesheet
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
        .findAll { e -> e.value }
        .collect { e -> e.value.toUriString().getChars() }
        .transpose()
        .findIndexOf { e ->
            def cs = e.unique()
            if (cs.size() != 1) return true
            c << cs.pop()
            return false
        }
    return file("${c.join('')}")
}
