#! /bin/python

## This utility script filters a multi-fasta file containing the consensus sequences of all protein families in a pangenome down to a multi-fasta file of all the core proteins using the pangenome presence/absence matrix.
##
## USAGE
## python extract_core_genes.py @core_threshold @pangenome_folder @output_folder
##
## PARAMS
## core_threshold   float between 0 and 1; the core genome threshold
## pangenome_folder path to the directory in which the presence/absence matrix and the pangenome multifasta are saved
## output_folder    path to the directory in which the extracted core protein sequences will be saved in a multi-fasta

import pandas as pd
from Bio import SeqIO
import sys

core_threshold = float(sys.argv[1])
pangenome_folder = sys.argv[2]
marker_folder = sys.argv[3]

matrix = pangenome_folder + "/matrix.csv"
prot_fams = pangenome_folder + "/all_protein_families.faa"

# Determine number of strains from the header of the matrix file
header = pd.read_table(matrix, sep = ",", nrows = 1)
n_strains = header.shape[1] - 14
n_strains_threshold = n_strains * core_threshold

# Only read the No. isolated column to get the strain abundance of the gene family
genes_pa_data = pd.read_table(matrix, sep = ",", usecols = [0,3])

# Read all pangenome protein sequences
with open(prot_fams, "r") as handle:
    genes_pan_ref_seqs = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    
# Filter
core_genes = genes_pa_data[genes_pa_data["No. isolates"] >= n_strains_threshold]["Gene"].to_list()
core_genes_seqs = [j for i,j in genes_pan_ref_seqs.items() if i in core_genes]

# Save filtered core protein sequences
with open(marker_folder + "/core_gene_seqs.faa", "w") as handle:
    SeqIO.write(core_genes_seqs, handle, "fasta")
