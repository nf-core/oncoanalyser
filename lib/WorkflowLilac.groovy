//
// This file holds several functions specific to the subworkflows/lilac.nf in the umccr/hmftools pipeline
//

class WorkflowLilac {

  public static get_slice_inputs(ch, hla_bed) {
    // channel: [val(meta_lilac), bam, bai, bed]
    def d = ch
      .flatMap { meta, tbam, nbam, tbai, nbai ->
        def sample_types = ['tumor': [tbam, tbai], 'normal': [nbam, nbai]]
        sample_types
          .keySet()
          .collect { sample_type ->
            def fps = sample_types[sample_type]
            def sample_name = meta.get(['sample_name', sample_type])
            def meta_lilac = [
              id: sample_name,
              sample_type: sample_type,
              meta_full: meta,
            ]
            return [meta_lilac, *fps, hla_bed]
          }
      }
    return d
  }

  public static get_unique_input_files(ch) {
    // channel: [val(meta_lilac), bam, bai, bed]
    def d = ch
      .map { [it[1..-1], it[0]] }
      .groupTuple()
      .map { filepaths, meta_lilacs ->
        def (meta_fulls, sample_types) = meta_lilacs
          .collect {
            [it.meta_full, it.sample_type]
          }
          .transpose()

        def sample_type = sample_types.unique(false)
        assert sample_type.size() == 1

        def id = meta_fulls.collect { it.id }.join('__')
        def meta_lilac_new = [
          id: "${id}_${sample_type[0]}",
          metas_full: meta_fulls,
          sample_type: sample_type[0],
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
        meta_lilac.metas_full.collect { meta -> [meta.id, meta, [meta_lilac.sample_type, *fps]] }
      }
      .groupTuple(size: 2)
      .map { id, meta, other ->
        def data = [:]
        other.each { sample_type, bam, bai ->
          data[[sample_type, 'bam']] = bam
          data[[sample_type, 'bai']] = bai
        }
        [
          meta[0],
          data.get(['tumor', 'bam']),
          data.get(['normal', 'bam']),
          data.get(['tumor', 'bai']),
          data.get(['normal', 'bai']),
        ]
      }
    return d
  }
}
