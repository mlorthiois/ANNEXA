process RESTORE_BIOTYPE {
  conda (params.enable_conda ? "conda-forge::python=3.10.4" : null)
  container 'quay.io/biocontainers/python:3.10.4'

  input:
  file ref
  file novel_isoforms

  output:
  path "novel.isoforms.gtf"

  """
  restore_ref_attributes.py -gtf ${novel_isoforms} -ref $ref > novel.isoforms.gtf
  """
}
