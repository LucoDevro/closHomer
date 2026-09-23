#! /bin/python

## This utility script adapts the network file produced by BiG-SCAPE and adds a classification column so that it can be easily imported into CytoScape.
##
## USAGE
## python annotate_bigscape_network.py @network_file @classification @level @antismash_output @out_file
##
## PARAMS
## network_file     the network file produced by BiG-SCAPE
## classification   classification table produced by the GTDB module
## level            classification level at which the grouping was defined
## antismash_output path to the folder with raw antiSMASH output
## out_file         location of the adapted network file

import pandas as pd
import os
import sys

bigscape_network_file = sys.argv[1]
classification_file = sys.argv[2]
level = sys.argv[3]
antismash_output = sys.argv[4]
output_file = sys.argv[5]

bigscape_network = pd.read_table(bigscape_network_file)
classification = pd.read_table(classification_file, usecols = ['accession', level], sep = "\t")

root = os.getcwd()
os.chdir(antismash_output)
group_annotations = {}

# Build the rRNA cluster metadata column by assigning the rRNA cluster to each detected antiSMASH region
for assembly in os.listdir():
    os.chdir(assembly)
    group = list(classification[classification['accession'] == assembly][level])[0]
    genbank_region_files = list(filter(lambda x: '.region0' in x, os.listdir()))
    for region in genbank_region_files:
        group_annotations['.'.join([assembly] + region.split('.')[:-1])] = group
    os.chdir('..')

# Assign an rRNA cluster to each node in the network file of BiG-SCAPE, for both sides of a link. BGCs that do not have a rRNA cluster annotation, usually are the MIBiG clusters.
bigscape_network_with_group = bigscape_network.copy(deep = True)
group1_annotations = []
for i in list(bigscape_network_with_group['Clustername 1']):
    try:
        group1_annotations.append(group_annotations[i])
    except KeyError:
        group1_annotations.append('MIBiG')
group2_annotations = []
for i in list(bigscape_network_with_group['Clustername 2']):
    try:
        group2_annotations.append(group_annotations[i])
    except KeyError:
        group2_annotations.append('MIBiG')
bigscape_network_with_group['group1'] = group1_annotations
bigscape_network_with_group['group2'] = group2_annotations

# Convert the unidirectional links into bidirectional ones by adding the reverse link
cols = list(bigscape_network_with_group.columns)
old_cols = cols.copy()
cols[0], cols[1], cols[13], cols[14] = cols[1], cols[0], cols[14], cols[13]
bigscape_network_with_group_swapped = bigscape_network_with_group[cols]
bigscape_network_with_group_swapped.columns = old_cols
bigscape_network_with_group_dupl = pd.concat([bigscape_network_with_group, bigscape_network_with_group_swapped]).reset_index(drop = True)

# Export new network file
os.chdir(root)
bigscape_network_with_group_dupl.to_csv(output_file, sep = "\t", index = False)
