process FILTER {
  conda (params.enable_conda ? "conda-forge::python=3.10.4" : null)
  container 'quay.io/biocontainers/python:3.10.4'
  publishDir "$params.outdir/final", mode: 'copy'

  input:
  file novel
  file tfkmers
  file bambu_ndr

  output:
  path "novel.filter.gtf"

  """
  filter_gtf_ndr.py \
    --gtf ${novel} \
    --tfkmers ${tfkmers} \
    --bambu ${bambu_ndr} \
    --tfkmers-threshold ${params.tfkmers_threshold} \
    --bambu-threshold ${params.bambu_threshold} \
    --operation ${params.operation} \
  | GTF.py format \
  > novel.filter.gtf
  """
}
