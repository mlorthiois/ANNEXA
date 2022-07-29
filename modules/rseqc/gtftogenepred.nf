process GTF_TO_GENEPRED {
  input:
  file gtf

  output:
  path "${gtf.simpleName}.genePred"

  """
  gtfToGenePred -genePredExt ${gtf} ${gtf.simpleName}.genePred
  """
}
