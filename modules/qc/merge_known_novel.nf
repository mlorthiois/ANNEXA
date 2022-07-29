process MERGE_ANNOTATIONS {
  publishDir "$params.outdir/final", mode: 'copy'

  input:
  file novel
  file ref
  val origin

  output:
  path "extended_annotations.${origin}.gtf"

  """
  cat ${novel} ${ref} | GTF.py format > extended_annotations.${origin}.gtf
  """
}
