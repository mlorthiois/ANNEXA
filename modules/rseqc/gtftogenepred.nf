process GTF_TO_GENEPRED {
  conda (params.enable_conda ? "bioconda::ucsc-gtftogenepred=377" : null)
  container 'quay.io/biocontainers/ucsc-gtftogenepred:377--ha8a8165_5'

  input:
  file gtf

  output:
  path "${gtf.simpleName}.genePred"

  """
  gtfToGenePred -genePredExt ${gtf} ${gtf.simpleName}.genePred
  """
}
