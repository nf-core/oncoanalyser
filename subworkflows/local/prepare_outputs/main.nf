//
// Publish helpers and output aggregation for the nf-core/oncoanalyser pipeline
//

nextflow.enable.types = true

include { getLongitudinalSampleName } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getNormalDnaSampleName    } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorDnaSampleName  } from '../utils_nfcore_oncoanalyser_pipeline/accessors'
include { getTumorRnaSampleName  } from '../utils_nfcore_oncoanalyser_pipeline/accessors'

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
        filepaths << ["${meta.case_id}/${dir}/${rel_path}", f]
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

workflow PREPARE_OUTPUTS {
    take:
    amber: Channel<Tuple<Map, Path>>
    bamtools_tumor: Channel<Tuple<Map, Path>>
    bamtools_normal: Channel<Tuple<Map, Path>>
    chord: Channel<Tuple<Map, Path>>
    cider: Channel<Tuple<Map, List<Path>>>
    cobalt: Channel<Tuple<Map, Path>>
    cuppa: Channel<Tuple<Map, Path>>
    esvee: Channel<Tuple<Map, Path>>
    align_rna_tumor: Channel<Tuple<Map, Path, Path>>
    isofox: Channel<Tuple<Map, Path>>
    lilac: Channel<Tuple<Map, Path>>
    linx_germline: Channel<Tuple<Map, Path>>
    linx_somatic: Channel<Tuple<Map, Path>>
    linx_somatic_visualiser: Channel<Tuple<Map, Path>>
    linxreport_html: Channel<Tuple<Map, Path>>
    multiqc: Channel<List<Path>>
    neo_annotated_fusions: Channel<Tuple<Map, Path>>
    neo_finder: Channel<Tuple<Map, Path>>
    neo_scorer_dir: Channel<Tuple<Map, Path>>
    orange_json: Channel<Tuple<Map, Path>>
    orange_pdf: Channel<Tuple<Map, Path>>
    pave_germline: Channel<Tuple<Map, Path>>
    pave_somatic: Channel<Tuple<Map, Path>>
    peach: Channel<Tuple<Map, Path>>
    purple: Channel<Tuple<Map, Path>>
    qsee: Channel<Tuple<Map, Path>>
    redux_tumor: Channel<Tuple<Map, Path>>
    redux_normal: Channel<Tuple<Map, Path>>
    redux_donor: Channel<Tuple<Map, Path>>
    sage_append_somatic: Channel<Tuple<Map, Path>>
    sage_append_germline: Channel<Tuple<Map, Path>>
    sage_germline: Channel<Tuple<Map, Path>>
    sage_somatic: Channel<Tuple<Map, Path>>
    sage_somatic_visualiser: Channel<Tuple<Map, Path>>
    sigs: Channel<Tuple<Map, Path>>
    teal_normal_bam: Channel<Tuple<Map, Path, Path>>
    teal_tumor_bam: Channel<Tuple<Map, Path, Path>>
    teal_tsvs: Channel<Tuple<Map, List<Path>>>
    virusbreakend_tsv: Channel<Tuple<Map, Path>>
    virusbreakend_vcf: Channel<Tuple<Map, Path>>
    virusinterpreter: Channel<Tuple<Map, Path>>
    write_reference_data: Channel<Path>
    command_files: Channel<Tuple<Map, String, List<Path>>>
    wisp: Channel<Tuple<Map, Path>>
    cobalt_normalisation_tsv: Channel<Path>
    isofox_normalisation_csv: Channel<Path>
    pave_pon_panel_creation_artefacts: Channel<Path>

    main:
    results = channel.empty()
        .mix(
            amber.flatMap { meta, d -> return get_dir_filepaths(meta, d) },
            bamtools_tumor.flatMap { meta, d -> return get_dir_filepaths(meta, d, "bamtools/${getTumorDnaSampleName(meta)}") },
            bamtools_normal.flatMap { meta, d -> return get_dir_filepaths(meta, d, "bamtools/${getNormalDnaSampleName(meta)}") },
            chord.flatMap { meta, d -> return get_dir_filepaths(meta, d) },
            cider.flatMap { meta, fps -> return fps.collect { d -> ["${meta.case_id}/cider/${d.name}", d] } },
            cobalt.flatMap { meta, d -> return get_dir_filepaths(meta, d) },
            cuppa.flatMap { meta, d -> return get_dir_filepaths(meta, d) },
            esvee.flatMap { meta, d -> return get_dir_filepaths(meta, d) },
            align_rna_tumor.flatMap { meta, bam, bai -> return [bam, bai].findAll().collect { d -> ["${meta.case_id}/alignments/${getTumorRnaSampleName(meta)}/${d.name}", d] } },
            isofox.flatMap { meta, d -> return get_dir_filepaths(meta, d) },
            lilac.flatMap { meta, d -> return get_dir_filepaths(meta, d) },
            linx_germline.flatMap { meta, d -> return get_dir_filepaths(meta, d, 'linx/germline_annotations') },
            linx_somatic.flatMap { meta, d -> return get_dir_filepaths(meta, d, 'linx/somatic_annotations') },
            linx_somatic_visualiser.flatMap { meta, d -> return get_dir_filepaths(meta, d, 'linx/somatic_plots') },
            linxreport_html.map { meta, d -> return ["${meta.case_id}/linx/${d.name}", d] },
            multiqc.flatMap { fps -> return fps.collect { d -> [d.name, d] } },
            neo_annotated_fusions.filter { meta, d -> d != null }.map { meta, d -> return ["${meta.case_id}/neo/annotated_fusions/${d.name}", d] },
            neo_finder.flatMap { meta, d -> return get_dir_filepaths(meta, d, 'neo/finder') },
            neo_scorer_dir.flatMap { meta, d -> return get_dir_filepaths(meta, d, 'neo/scorer') },
            orange_json.filter { meta, d -> d != null }.map { meta, d -> return ["${meta.case_id}/orange/${d.name}", d] },
            orange_pdf.filter { meta, d -> d != null }.map { meta, d -> return ["${meta.case_id}/orange/${d.name}", d] },
            pave_germline.flatMap { meta, d -> return get_dir_filepaths(meta, d, 'pave/germline') },
            pave_somatic.flatMap { meta, d -> return get_dir_filepaths(meta, d, 'pave/somatic') },
            peach.flatMap { meta, d -> return get_dir_filepaths(meta, d) },
            purple.flatMap { meta, d -> return get_dir_filepaths(meta, d) },
            qsee.flatMap { meta, d -> return get_dir_filepaths(meta, d) },
            redux_tumor.flatMap { meta, d -> return get_dir_filepaths(meta, d, "alignments/${getTumorDnaSampleName(meta)}") },
            redux_normal.flatMap { meta, d -> return get_dir_filepaths(meta, d, "alignments/${getNormalDnaSampleName(meta)}") },
            redux_donor.flatMap { meta, dirs ->
                def donor_dirs = dirs == null ? [] : (dirs instanceof List ? dirs : [dirs])
                return donor_dirs.collectMany { d -> get_dir_filepaths(meta, d, "alignments/${d.name.replaceFirst(/^redux_/, '')}") }
            },
            sage_append_somatic.flatMap { meta, d -> return get_dir_filepaths(meta, d, "sage_append/${getLongitudinalSampleName(meta) ?: getTumorDnaSampleName(meta)}") },
            sage_append_germline.flatMap { meta, d -> return get_dir_filepaths(meta, d, "sage_append/${getNormalDnaSampleName(meta)}") },
            sage_germline.flatMap { meta, d -> return get_dir_filepaths(meta, d, 'sage/germline') },
            sage_somatic.flatMap { meta, d -> return get_dir_filepaths(meta, d, 'sage/somatic') },
            sage_somatic_visualiser.flatMap { meta, d -> return get_dir_filepaths(meta, d, 'sage/visualiser') },
            sigs.flatMap { meta, d -> return get_dir_filepaths(meta, d) },
            teal_normal_bam.flatMap { meta, bam, bai -> return [bam, bai].findAll().collect { d -> ["${meta.case_id}/teal/${d.name}", d] } },
            teal_tumor_bam.flatMap { meta, bam, bai -> return [bam, bai].findAll().collect { d -> ["${meta.case_id}/teal/${d.name}", d] } },
            teal_tsvs.flatMap { meta, e -> def fps = e instanceof Collection ? e : [e]; fps.collect { d -> ["${meta.case_id}/teal/${d.name}", d] } },
            virusbreakend_tsv.map { meta, d -> return ["${meta.case_id}/virusbreakend/${d.name}", d] },
            virusbreakend_vcf.map { meta, d -> return ["${meta.case_id}/virusbreakend/${d.name}", d] },
            virusinterpreter.flatMap { meta, d -> return get_dir_filepaths(meta, d) },
            wisp.flatMap { meta, d -> return get_dir_filepaths(meta, d) },
            write_reference_data.map { d -> return ["reference_data/${workflow.manifest.version}/", d] },
            cobalt_normalisation_tsv.map { d -> return ["panel_resources/${d.name}", d] },
            isofox_normalisation_csv.map { d -> return ["panel_resources/${d.name}", d] },
            pave_pon_panel_creation_artefacts.map { d -> return ["panel_resources/${d.name}", d] },
            command_files.flatMap { f -> get_command_log_filepath(f) }
        )
        .flatMap { meta, d -> return d instanceof Collection ? d.collect { e -> [meta, e] } : [[meta, d]] }

    emit:
    results = results
}
