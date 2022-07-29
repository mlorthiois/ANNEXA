process GENEPRED_TO_BED {
  input:
  file genePred

  output:
  path "${genePred.simpleName}.bed12"

  """
  genePredToBed ${genePred} ${genePred.simpleName}.bed12
  """
}
