process INDEX_BAM {
  input:
  file bam

  output:
  file("*.bai")

  """
  samtools index $bam
  """
}
