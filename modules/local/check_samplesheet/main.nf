process CHECK_SAMPLESHEET {
  conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
  container 'quay.io/biocontainers/python:3.9--1'

  input:
  path samplesheet
  val mode

  output:
  path samplesheet

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''

  """
  check_samplesheet.py \\
    --input_fp "${samplesheet}" \\
    --mode "${mode}"

  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      python: \$(python --version | sed 's/Python //g')
  END_VERSIONS
  """

  stub:
  """
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
