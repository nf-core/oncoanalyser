//
// This file holds several functions specific to the subworkflows/lilac.nf in the nf-core/oncoanalyser pipeline
//
import Constants
import Utils


class WorkflowLilac {

    public static getSliceInputs(ch, ch_slice_bed) {
        // channel: [val(meta_lilac), bam, bai]
        def d = ch
            .flatMap { meta, tbam, nbam, tbai, nbai ->
                def sample_types = [Constants.DataType.TUMOR, Constants.DataType.NORMAL]
                sample_types
                    .collect { sample_type ->
                        def fps
                        if (sample_type == Constants.DataType.TUMOR) {
                            fps = [tbam, tbai]
                        } else if (sample_type == Constants.DataType.NORMAL) {
                            fps = [nbam, nbai]
                        } else {
                            assert false : "got bad sample type"
                        }
                        def sample_name = meta.get(['sample_name', sample_type])
                        def meta_lilac = [
                            key: meta.id,
                            id: sample_name,
                            // NOTE(SW): must use string representation for caching purposes
                            sample_type_str: sample_type.name(),
                        ]
                        return [meta_lilac, *fps]
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
                def (keys, sample_names, sample_type_strs) = meta_lilac
                    .collect {
                        [it.key, it.id, it.sample_type_str]
                    }
                    .transpose()

                def sample_type_strs_unique = sample_type_strs.unique(false)
                assert sample_type_strs_unique.size() == 1
                def sample_type_str = sample_type_strs_unique[0]

                def key = keys.join('__')
                def meta_lilac_new = [
                    keys: keys,
                    id: sample_names.join('__'),
                    id_simple: keys.join('__'),
                    sample_type_str: sample_type_str,
                ]
                return [meta_lilac_new, *filepaths]
            }
        return d
    }

    public static sortSlices(ch) {
        // Collect T/N pairs into single channel element
        // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai]
        def d = ch
            .flatMap{ data ->
                def meta_lilac = data[0]
                def fps = data[1..-1]
                meta_lilac.keys.collect { key ->
                    return [key, [meta_lilac.sample_type_str, *fps]]
                }
            }
            .groupTuple(size: 2)
            .map { key, other ->
                def data = [:]
                other.each { sample_type_str, bam, bai ->
                    def sample_type = Utils.getEnumFromString(sample_type_str, Constants.DataType)
                    data[[sample_type, 'bam']] = bam
                    data[[sample_type, 'bai']] = bai
                }
                [
                    [key: key],
                    data.get([Constants.DataType.TUMOR, 'bam']),
                    data.get([Constants.DataType.NORMAL, 'bam']),
                    data.get([Constants.DataType.TUMOR, 'bai']),
                    data.get([Constants.DataType.NORMAL, 'bai']),
                ]
            }
        return d
    }
}
