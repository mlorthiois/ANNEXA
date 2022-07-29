process VALIDATE_INPUT_GTF {
  input:
  file gtf

  output:
  file 'input.formatted.gtf'

  """
  validate_gtf.py ${gtf} > input.formatted.gtf
  """
}
