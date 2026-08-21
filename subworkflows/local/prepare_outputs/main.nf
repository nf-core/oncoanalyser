//
// Publish helpers shared by the workflow emit blocks
//

// NOTE(LN): Nextflow bug as of NXF_VER=26.04.6
//
// When using google cloud executor (and potentially other cloud executors):
// - Process emits output as directory (e.g. amber_dir) -> directory publish FAIL -> BUG
// - Process emits output as files or glob (e.g. teal_tsvs) -> individual files publish SUCCESS
//
// To publish directories, the workaround is to recursively extract the file paths within that directory
//
def get_dir_filepaths(meta, d, target_dir=null) {

    // NOTE(SW): restored channels carry [meta, null] for skipped cases; a topic would just not emit
    if (d == null) { return [] }

    def dir = target_dir ?: d.name

    def filepaths = []
    d.eachFileRecurse(groovy.io.FileType.FILES) { f ->
        def rel_path = d.relativize(f).toString()
        filepaths << ["${meta.key ?: meta.case_id}/${dir}/${rel_path}", f]
    }
    return filepaths
}


def get_command_log_filepath(data) {

    def other_logs = ['gatk4_bwa_index_image', 'gridss_index', 'bwa_index', 'bwamem2_index', 'samtools_dict', 'samtools_faidx', 'star_genomegenerate', 'multiqc']
    def panel_logs = ['cobalt_panel_normalisation', 'pave_pon_panel_creation', 'isofox_panel_normalisation']

    def (meta, name, fps_all) = data

    def fps = fps_all.findAll { f -> f.name.matches(/.*\.command\.(sh|out|err|log|run)/) }

    if (name == 'extracttarball') {
        return fps.collect { d -> ["logs/other/${name}.${meta.id}${d.name}", d] }
    } else if (other_logs.contains(name)) {
        return fps.collect { d -> ["logs/other/${name}${d.name}", d] }
    } else if (panel_logs.contains(name)) {
        return fps.collect { d -> ["logs/panel_resources/${name}${d.name}", d] }
    } else {
        return fps.collect { d -> ["logs/${meta.key}/${name}.${meta.id}${d.name}", d] }
    }

}
