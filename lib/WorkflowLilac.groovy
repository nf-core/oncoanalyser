//
// This file holds several functions specific to the subworkflows/lilac.nf in the nf-core/oncoanalyser pipeline
//

class WorkflowLilac {

    public static get_slice_inputs(ch, ch_slice_bed) {
        // channel: [val(meta_lilac), bam, bai]
        def d = ch
            .flatMap { meta, tbam, nbam, tbai, nbai ->
                def sample_types = ['tumor': [tbam, tbai], 'normal': [nbam, nbai]]
                sample_types
                    .keySet()
                    .collect { sample_type ->
                        def fps = sample_types[sample_type]
                        def sample_name = meta.get(['sample_name', sample_type])
                        def meta_lilac = [
                            key: meta.id,
                            id: sample_name,
                            sample_type: sample_type,
                        ]
                        return [meta_lilac, *fps]
                    }
            }
        // channel: [val(meta_lilac), bam, bai, bed]
        return d.combine(ch_slice_bed)
    }

    public static get_unique_input_files(ch) {
        // channel: [val(meta_lilac), bam, bai, bed]
        def d = ch
            .map { [it[1..-1], it[0]] }
            .groupTuple()
            .map { filepaths, meta_lilac ->
                def (keys, sample_names, sample_types) = meta_lilac
                    .collect {
                        [it.key, it.id, it.sample_type]
                    }
                    .transpose()

                def sample_types_unique = sample_types.unique(false)
                assert sample_types_unique.size() == 1
                def sample_type = sample_types_unique[0]

                def key = keys.join('__')
                def meta_lilac_new = [
                    keys: keys,
                    id: sample_names.join('__'),
                    id_simple: keys.join('__'),
                    sample_type: sample_type,
                ]
                return [meta_lilac_new, *filepaths]
            }
        return d
    }

    public static sort_slices(ch) {
        // Collect T/N pairs into single channel element
        // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai]
        def d = ch
            .flatMap{ data ->
                def meta_lilac = data[0]
                def fps = data[1..-1]
                meta_lilac.keys.collect { key ->
                    return [key, [meta_lilac.sample_type, *fps]]
                }
            }
            .groupTuple(size: 2)
            .map { key, other ->
                def data = [:]
                other.each { sample_type, bam, bai ->
                    data[[sample_type, 'bam']] = bam
                    data[[sample_type, 'bai']] = bai
                }
                [
                    [key: key],
                    data.get(['tumor', 'bam']),
                    data.get(['normal', 'bam']),
                    data.get(['tumor', 'bai']),
                    data.get(['normal', 'bai']),
                ]
            }
        return d
    }
}
