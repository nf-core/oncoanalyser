process VIRUSINTERPRETER {
  container 'docker.io/scwatts/virus_interpreter:1.2--0'

  input:
  tuple val(meta), path(virus_tsv), path(purple_purity), path(purple_qc), path(wgs_metrics)
  path taxonomy
  path virus_reporting

  output:
  tuple val(meta), path('virusinterpreter'), emit: virusinterpreter_dir
  path 'versions.yml'                      , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''

  """
  mkdir -p virusinterpreter/

  java \\
    -Xmx${task.memory.giga}g \\
    -jar ${task.ext.jarPath} \\
      -sample_id ${meta.get(['sample_name', 'tumor'])} \\
      -purple_purity_tsv ${purple_purity} \\
      -purple_qc_file ${purple_qc} \\
      -tumor_sample_wgs_metrics_file ${wgs_metrics} \\
      -virus_breakend_tsv ${virus_tsv} \\
      -taxonomy_db_tsv ${taxonomy} \\
      -virus_reporting_db_tsv ${virus_reporting} \\
      -output_dir ./virusinterpreter/

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      "virus interpreter": \$(java -jar "${task.ext.jarPath}" | sed -n '1s/^.*Interpreter v//p')
  END_VERSIONS
  """

  stub:
  """
  mkdir virusinterpreter/
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
