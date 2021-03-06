#! /usr/bin/env python3
import argparse
from GTF import GTF

#######################################################
parser = argparse.ArgumentParser(description="Utility tools for your GTF files.")
parser.add_argument(
    "-mRNA", help="Path to your mRNA only GTF file.", type=argparse.FileType("r"),
)
parser.add_argument(
    "-lncRNA", help="Path to your lncRNA only GTF file.", type=argparse.FileType("r"),
)
args = parser.parse_args()

#######################################################
# Parse gtf classed as protein_coding by FEELnc
gene_ids = {}
# Add transcript_biotype = protein_coding
for record in GTF.parse_by_line(args.mRNA):
    record.source = "FEELnc"
    record["transcript_biotype"] = "protein_coding"
    record["gene_biotype"] = "protein_coding"

    if record["gene_id"] in gene_ids:
        gene_ids[record["gene_id"]].append(record)
    else:
        gene_ids[record["gene_id"]] = [record]

#######################################################
# Parse gtf classed as lncRNA by FEELnc
lncRNAs = {}
# If gene_id also found in protein_coding => protein_coding else lncRNA
for record in GTF.parse_by_line(args.lncRNA):
    record.source = "FEELnc"
    if record["gene_id"] in gene_ids:
        record["transcript_biotype"] = "protein_coding"
        record["gene_biotype"] = "protein_coding"
        gene_ids[record["gene_id"]].append(record)

    else:
        record["transcript_biotype"] = "lncRNA"
        record["gene_biotype"] = "lncRNA"

        if record["gene_id"] in lncRNAs:
            lncRNAs[record["gene_id"]].append(record)
        else:
            lncRNAs[record["gene_id"]] = [record]

#######################################################
# Keep only these attributes
filter = [
    "gene_id",
    "transcript_id",
    "exon_number",
    "gene_biotype",
    "transcript_biotype",
]

# Print to stdout gtf
for gene in gene_ids:
    for exon in gene_ids[gene]:
        exon.filter_attributes(filter)
        print(exon)

for gene in lncRNAs:
    for exon in lncRNAs[gene]:
        exon.filter_attributes(filter)
        print(exon)
