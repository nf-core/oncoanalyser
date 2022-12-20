//
// This file holds several functions specific to the subworkflows/lilac.nf in the nf-core/oncoanalyser pipeline
//
import nextflow.Nextflow

import Constants
import Utils


class WorkflowLilac {

    public static getSliceInputs(ch, ch_slice_bed) {
        // channel: [val(meta_lilac), bam, bai]
        def d = ch
            .flatMap { meta, nbam_wgs, nbai_wgs, tbam_wgs, tbai_wgs, tbam_wts, tbai_wts ->
                def data = [
                    [nbam_wgs, nbai_wgs, Constants.FileType.BAM_WGS, Constants.DataType.NORMAL],
                    [tbam_wgs, tbai_wgs, Constants.FileType.BAM_WGS, Constants.DataType.TUMOR],
                    [tbam_wts, tbai_wts, Constants.FileType.BAM_WTS, Constants.DataType.TUMOR],
                ]
                def data_present = data.findAll { it[0] }
                data_present
                    .collect { bam, bai, filetype, sample_type ->
                        def sample_name = meta.get(['sample_name', sample_type])
                        def meta_lilac = [
                            key: meta.id,
                            id: sample_name,
                            // NOTE(SW): must use string representation for caching purposes
                            sample_type_str: sample_type.name(),
                            filetype_str: filetype.name(),
                            count: data_present.size(),
                        ]
                        return [meta_lilac, bam, bai]
                    }
            }
        // channel: [val(meta_lilac), bam, bai, bed]
        return d.combine(ch_slice_bed)
    }

    public static getUniqueInputFiles(ch) {
        // channel: [val(meta_lilac), bam, bai, bed]
        def d = ch
            .map { [it[1..-1], it[0]] }
            .groupTuple()
            .map { filepaths, meta_lilac ->
                // NOTE(SW): pattern needs to be generalised
                def (keys, sample_names, filetype_strs, sample_type_strs) = meta_lilac
                    .collect { [it.key, it.id, it.filetype_str, it.sample_type_str] }
                    .transpose()
                def sample_type_str = getValue(sample_type_strs)
                def filetype_str = getValue(filetype_strs)

                def key = keys.join('__')
                def meta_lilac_new = [
                    keys: keys,
                    id: sample_names.join('__'),
                    id_simple: keys.join('__'),
                    filetype_str: filetype_str,
                    sample_type_str: sample_type_str,
                ]
                return [meta_lilac_new, *filepaths]
            }
        return d
    }

    public static sortSlices(ch) {
        // Gather and order files from same grouping using non-blocked groupTuple via provided file counts
        // channel: [val(meta_lilac), normal_wgs_bam, normal_wgs_bai, tumor_wgs_bam, tumor_wgs_bai, tumor_wts_bam, tumor_wts_bai]
        def d = ch
            .map { meta, data ->
                return [nextflow.Nextflow.groupKey(meta.key, meta.count), meta, data]
            }
            .groupTuple()
            .map { gk, metas, values ->
                assert metas.unique().size() == 1
                def meta = metas[0]
                def data = [:]
                values.each { filetype_str, sample_type_str, bam, bai ->
                    def sample_type = Utils.getEnumFromString(sample_type_str, Constants.DataType)
                    def filetype = Utils.getEnumFromString(filetype_str, Constants.FileType)
                    data[[sample_type, filetype, 'bam']] = bam
                    data[[sample_type, filetype, 'bai']] = bai
                }
                return [
                    meta,
                    data.get([Constants.DataType.NORMAL, Constants.FileType.BAM_WGS, 'bam'], []),
                    data.get([Constants.DataType.NORMAL, Constants.FileType.BAM_WGS, 'bai'], []),
                    data.get([Constants.DataType.TUMOR, Constants.FileType.BAM_WGS, 'bam'], []),
                    data.get([Constants.DataType.TUMOR, Constants.FileType.BAM_WGS, 'bai'], []),
                    data.get([Constants.DataType.TUMOR, Constants.FileType.BAM_WTS, 'bam'], []),
                    data.get([Constants.DataType.TUMOR, Constants.FileType.BAM_WTS, 'bai'], []),
                ]
            }
        return d
    }

    public static getValue(l) {
        def u = l.unique(false)
        assert u.size() == 1
        return u[0]
    }
}
