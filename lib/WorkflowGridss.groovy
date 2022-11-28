//
// This file holds several functions specific to the subworkflows/gridss.nf in the umccr/hmftools pipeline
//

class WorkflowGridss {

  public static get_inputs(ch) {
    // channel (a): [val(meta_gridss), bam, bai, sv]
    // channel (b): [val(meta_gridss), bam]
    def sample_types = ['tumor', 'normal']
    def d = ch
      .flatMap { meta ->
        sample_types
          .collect { sample_type ->
            def key_sname = ['sample_name', sample_type]
            def key_bam = ['bam_wgs', sample_type]
            def key_sv = ['vcf_sv', sample_type]
            if (! meta.containsKey(key_sname) && ! meta.containsKey(key_bam)) {
              return []
            }
            def meta_gridss = [
              id: meta.get(key_sname),
              group_key: meta.id,
              sample_type: sample_type,
              subject_name: meta.subject_name,
            ]
            def v = [meta_gridss, meta.get(key_bam)]
            def has_sv = meta.containsKey(key_sv)
            if (has_sv) {
              v = v + ["${meta.get(key_bam)}.bai", meta.get(key_sv)]
            }
            return [has_sv, v]
          }
      }
      return d
  }

  public static get_unique_input_files(ch) {
    // channel (a): [val(meta_gridss), bam, bai, sv]
    // channel (b): [val(meta_gridss), bam]
    def d = ch
      .map { [it[1..-1], it[0]] }
      // NOTE(SW): number of grouped elements is unknown here but does not block since all inputs
      // are derived directly from the samplesheet
      .groupTuple()
      .map { filepaths, meta_gridsss ->
        def (sample_names, ids, sample_types, subject_names) = meta_gridsss
          .collect {
            [it.id, it.group_key, it.sample_type, it.subject_name]
          }
          .transpose()

        def sample_name = sample_names.unique(false)
        def sample_type = sample_types.unique(false)
        def subject_name = subject_names.unique(false)
        assert sample_name.size() == 1
        assert sample_type.size() == 1
        assert subject_name.size() == 1

        def meta_gridss_new = [
          id: ids.join('__'),
          id_list: ids,
          sample_name: sample_name[0],
          sample_type: sample_type[0],
          subject_name: subject_name[0],
        ]
        return [meta_gridss_new, *filepaths]
      }
    return d
  }

  public static get_assemble_inputs(ch) {
    // NOTE(SW): elements in the bams, preprocess_dirs, and labels are linked with via consistent ordering
    // channel: [val(meta_gridss), [bams], [preprocess_dirs], [labels]]
    def d = ch
      .map { subject_name, other ->
        def data = [:]
        def ids = [] as Set

        // NOTE(SW): must determine whether input ordering is determininistic here

        other
          .each { meta_gridss, bam, preprocess_dir ->
            data.get([meta_gridss.sample_type, 'label'], []) << meta_gridss.sample_name
            data.get([meta_gridss.sample_type, 'bam_wgs'], []) << bam
            data.get([meta_gridss.sample_type, 'preprocess_dir'], []) << preprocess_dir
            ids += meta_gridss.id_list
          }
        def meta_gridss_new = [
          id: ids.join('__'),
          id_list: ids,
        ]

        def sample_types = ['normal', 'tumor']
        return [
          meta_gridss_new,
          sample_types.inject([]) { d, sample_type -> d + data.get([sample_type, 'bam_wgs'], []) },
          sample_types.inject([]) { d, sample_type -> d + data.get([sample_type, 'preprocess_dir'], []) },
          sample_types.inject([]) { d, sample_type -> d + data.get([sample_type, 'label'], []) },
        ]
      }
    return d
  }
}
