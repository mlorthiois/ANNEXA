process FILTER {
  publishDir "$params.outdir/final", mode: 'copy'

  input:
  file novel
  file mouloud
  file bambu_ndr

  output:
  path "novel.filter.gtf"

  """
  mv novel.gtf novel.filter.gtf
  """
}
