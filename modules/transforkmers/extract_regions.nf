process EXTRACT_TSS_REGIONS {
  input:
  file novel_gtf

  output:
  path "tss_regions.bed6", emit: tss_regions

  """
  cat ${novel_gtf} | extract_tss.py -l 512 > tss_regions.bed6
  """
}
