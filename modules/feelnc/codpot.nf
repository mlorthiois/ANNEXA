process FEELNC_CODPOT {
  memory params.maxMemory

  input:
  file ref
  file fa
  file gtf

  output:
  path "feelnc_codpot_out/new.lncRNA.gtf", emit: lncRNA
  path "feelnc_codpot_out/new.mRNA.gtf", emit: mRNA

  """
  grep "protein_coding" ${ref} > known_mRNA.gtf
  grep -v "protein_coding" ${ref} > known_lncRNA.gtf
  FEELnc_codpot.pl \
      -i $gtf \
      -g $fa \
      -a known_mRNA.gtf \
      -l known_lncRNA.gtf \
      --numtx=3000,3000 \
      -o new
  """
}

