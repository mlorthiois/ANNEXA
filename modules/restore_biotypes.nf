process RESTORE_BIOTYPE {
  input:
  file ref
  file novel_isoforms

  output:
  path "known.restored.gtf"

  """
  restore_ref_attributes.py -gtf ${novel_isoforms} -ref $ref > known.restored.gtf
  """
}
