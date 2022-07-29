process EXTRACT_TSS_SEQUENCES {
  input:
  file tss_regions
  file fa

  output:
  path "tss.fa", emit: tss_sequences

  """
  bedtools getfasta -nameOnly -fi ${fa} -bed ${tss_regions} | tr a-z A-Z > tss.fa
  """
}
