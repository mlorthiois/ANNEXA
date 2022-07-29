process MERGE_NOVEL {
  input:
  file novel_genes
  file novel_isoforms

  output:
  path "novel.gtf"

  """
  cat ${novel_genes} ${novel_isoforms} | GTF.py format > novel.gtf
  """
}
