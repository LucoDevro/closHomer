## This utility script replaces the RefSeq Assembly accession identifiers in the phylogeny by their species name, as linked in the metadata file
##
## USAGE
## python annotateContrees.py @tree @metadata @filename_replacement
##
## PARAMETERS
## tree                     a phyloheny in Newick format, with RefSeq Assembly accession identifiers as tip labels
## metadata                 metadata table in tsv format, containing the link between the identifiers (column "Genome_accession") and species names (column "Species")
## filename_replacement     part of the filename to append a "_cleaned" suffix to distinguish it from the original file, usually simply the file extension (e.g. ".contree")

import pandas as pd
import sys

tree_filename = sys.argv[1]
meta_filename = sys.argv[2]
str_repl = sys.argv[3]
names_by_ids = pd.read_table(meta_filename).set_index('Genome_accession').to_dict()['Species']

with open(tree_filename,'r') as handle:
    tree = handle.read()
for strain in names_by_ids.keys():
    tree = tree.replace(strain, names_by_ids[strain])
with open(tree_filename.replace(str_repl, str_repl + '_cleaned'), 'w') as handle:
    handle.write(tree)
