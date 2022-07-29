process FEELNC_FORMAT {
  publishDir "$params.outdir/final", mode: 'copy'

  input:
  file mRNA
  file lncRNA

  output:
  path "novel.full.gtf"

  """
  merge_feelnc.py --lncRNA ${lncRNA} --mRNA ${mRNA} > novel.full.gtf
  """
}
