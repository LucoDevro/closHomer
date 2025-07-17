## This utility script replaces the RefSeq Assembly accession identifiers in the phylogeny by their species name, as linked in the metadata file
##
## USAGE
## python annotate_tree.py @tree @metadata @filename_replacement
##
## PARAMETERS
## tree                     a phyloheny in Newick format, with RefSeq Assembly accession identifiers as tip labels
## metadata                 metadata table in tsv format, containing the link between the identifiers (column "Genome_accession") and species names (column "Species")
## filename_replacement     part of the filename to append a "_cleaned" suffix to distinguish it from the original file, usually simply the file extension (e.g. ".contree")

import pandas as pd
import sys
import re

tree_filename = sys.argv[1]
meta_filename = sys.argv[2]
str_repl = sys.argv[3]

def clean_species_name(name):
    new_name = re.sub(r"[\[\]]", "", name)
    new_name = re.sub(r"\(.*\)", "", new_name)
    new_name = re.sub(r" = .*", "", new_name)
    return new_name.rstrip()

metadata = pd.read_table(meta_filename, converters={'Organism Name': clean_species_name}).fillna('').set_index('Assembly Accession')
metadata['Full Organism Name'] = metadata['Organism Name'] + ' ' + metadata['Organism Infraspecific Names Strain']
metadata['Full Organism Name'] = metadata['Full Organism Name'].str.rstrip()
duplicate_species = metadata[metadata.duplicated(subset = "Organism Name", keep = False)]
metadata.loc[duplicate_species.index, "Organism Name"] = duplicate_species['Full Organism Name']
names_by_ids = metadata.to_dict()['Full Organism Name']

with open(tree_filename,'r') as handle:
    tree = handle.read()
for strain in names_by_ids.keys():
    tree = tree.replace(strain, names_by_ids[strain])
with open(tree_filename.replace(str_repl, str_repl + '_cleaned'), 'w') as handle:
    handle.write(tree)
