process FEELNC_FORMAT {
  conda (params.enable_conda ? "conda-forge::python=3.10.4" : null)
  container 'quay.io/biocontainers/python:3.10.4'

  input:
  file mRNA
  file lncRNA

  output:
  path "novel.genes.gtf"

  """
  merge_feelnc.py --lncRNA ${lncRNA} --mRNA ${mRNA} > novel.genes.gtf
  """
}
