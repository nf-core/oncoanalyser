//
// This file holds several functions specific to the subworkflows/lilac.nf in the nf-core/oncoanalyser pipeline
//
import nextflow.Nextflow

import Constants
import Utils


class WorkflowLilac {

    public static getSliceInputs(ch, tumor_sequence_type) {
        // channel: [val(meta_lilac), bam, bai]
        return ch
            .flatMap { meta, nbam_wgs, nbai_wgs, tbam_wgs, tbai_wgs, tbam_wts, tbai_wts ->
                def data = [
                    [nbam_wgs, nbai_wgs, Constants.SequenceType.WGS, Constants.SampleType.NORMAL],
                    [tbam_wgs, tbai_wgs, tumor_sequence_type, Constants.SampleType.TUMOR],
                    [tbam_wts, tbai_wts, Constants.SequenceType.WTS, Constants.SampleType.TUMOR],
                ]
                def data_present = data.findAll { it[0] }
                data_present
                    .collect { bam, bai, sequence_type, sample_type ->
                        def sample_name = meta.getAt(['sample_name', sample_type, sequence_type])
                        def meta_lilac = [
                            key: meta.id,
                            id: sample_name,
                            // NOTE(SW): must use string representation for caching purposes
                            sample_type_str: sample_type.name(),
                            sequence_type_str: sequence_type.name(),
                            count: data_present.size(),
                        ]
                        return [meta_lilac, bam, bai]
                    }
            }
    }

    public static getUniqueInputFiles(ch) {
        // channel: [val(meta_lilac), bam, bai]
        def d = ch
            .map { [it[1..-1], it[0]] }
            .groupTuple()
            .map { filepaths, meta_lilac ->
                // NOTE(SW): pattern needs to be generalised
                def (keys, sample_names, sequence_type_strs, sample_type_strs) = meta_lilac
                    .collect { [it.key, it.id, it.sequence_type_str, it.sample_type_str] }
                    .transpose()
                def sample_type_str = getValue(sample_type_strs)
                def sequence_type_str = getValue(sequence_type_strs)

                def key = keys.join('__')
                def meta_lilac_new = [
                    keys: keys,
                    id: sample_names.join('__'),
                    id_simple: keys.join('__'),
                    sequence_type_str: sequence_type_str,
                    sample_type_str: sample_type_str,
                ]
                return [meta_lilac_new, *filepaths]
            }
        return d
    }

    public static sortSlices(ch, tumor_sequence_type) {
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
                values.each { sequence_type_str, sample_type_str, bam, bai ->
                    def sample_type = Utils.getEnumFromString(sample_type_str, Constants.SampleType)
                    def sequence_type = Utils.getEnumFromString(sequence_type_str, Constants.SequenceType)
                    data[[sample_type, sequence_type, 'bam']] = bam
                    data[[sample_type, sequence_type, 'bai']] = bai
                }
                return [
                    meta,
                    data.get([Constants.SampleType.NORMAL, Constants.SequenceType.WGS, 'bam'], []),
                    data.get([Constants.SampleType.NORMAL, Constants.SequenceType.WGS, 'bai'], []),
                    data.get([Constants.SampleType.TUMOR, tumor_sequence_type, 'bam'], []),
                    data.get([Constants.SampleType.TUMOR, tumor_sequence_type, 'bai'], []),
                    data.get([Constants.SampleType.TUMOR, Constants.SequenceType.WTS, 'bam'], []),
                    data.get([Constants.SampleType.TUMOR, Constants.SequenceType.WTS, 'bai'], []),
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
