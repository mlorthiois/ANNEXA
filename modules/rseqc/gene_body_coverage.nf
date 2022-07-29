process GENEBODY_COVERAGE {
  publishDir "${params.outdir}/qc/${prefix}/rseqc", mode: 'copy'

  input:
  file bed 
  file "*"
  file "*"
  val prefix

  output:
  file "${prefix}.${bed.simpleName}.geneBodyCoverage.curves.pdf"
  file "${prefix}.${bed.simpleName}.geneBodyCoverage.txt"

  """
  geneBody_coverage.py \
    -i ./ \
    -r ${bed} \
    -o ${prefix}.${bed.simpleName}
  """
}
