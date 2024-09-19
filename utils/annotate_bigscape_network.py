#! /bin/python

## This utility script adapts the network file produced by BiG-SCAPE and adds an rRNA cluster annotation column so that it can be easily imported into CytoScape.
##
## USAGE
## python annotate_bigscape_network.py @network_file @metadata @antismash_output @out_file
##
## PARAMS
## network_file     the network file produced by BiG-SCAPE
## metadata         metadata of the original genome set
## antismash_output path to the folder with raw antiSMASH output
## out_file         location of the adapted network file

import pandas as pd
import os
import sys

bigscape_network_file = sys.argv[1]
metadata_file = sys.argv[2]
antismash_output = sys.argv[3]
output_file = sys.argv[4]

bigscape_network = pd.read_table(bigscape_network_file)
metadata = pd.read_table(metadata_file)

root = os.getcwd()
os.chdir(antismash_output)
cluster_annotations = {}

# Build the rRNA cluster metadata column by assigning the rRNA cluster to each detected antiSMASH region
for assembly in os.listdir():
    os.chdir(assembly)
    cluster = list(metadata[metadata['Genome_accession'] == assembly]['Taxonomic_cluster'])[0]
    genbank_region_files = list(filter(lambda x: '.region0' in x, os.listdir()))
    for region in genbank_region_files:
        cluster_annotations['.'.join(region.split('.')[:-1])] = cluster
    os.chdir('..')

# Assign an rRNA cluster to each node in the network file of BiG-SCAPE, for both sides of a link. BGCs that do not have a rRNA cluster annotation, usually are the MIBiG clusters.
bigscape_network_with_cluster = bigscape_network.copy(deep = True)
cluster1_annotations = []
for i in list(bigscape_network_with_cluster['Clustername 1']):
    try:
        cluster1_annotations.append(cluster_annotations[i])
    except KeyError:
        cluster1_annotations.append('MIBiG')
cluster2_annotations = []
for i in list(bigscape_network_with_cluster['Clustername 2']):
    try:
        cluster2_annotations.append(cluster_annotations[i])
    except KeyError:
        cluster2_annotations.append('MIBiG')
bigscape_network_with_cluster['cluster1'] = cluster1_annotations
bigscape_network_with_cluster['cluster2'] = cluster2_annotations

# Convert the unidirectional links into bidirectional ones by adding the reverse link
cols = list(bigscape_network_with_cluster.columns)
old_cols = cols.copy()
cols[0], cols[1], cols[13], cols[14] = cols[1], cols[0], cols[14], cols[13]
bigscape_network_with_cluster_swapped = bigscape_network_with_cluster[cols]
bigscape_network_with_cluster_swapped.columns = old_cols
bigscape_network_with_cluster_dupl = pd.concat([bigscape_network_with_cluster, bigscape_network_with_cluster_swapped]).reset_index(drop = True)

# Export new network file
os.chdir(root)
bigscape_network_with_cluster_dupl.to_csv(output_file, sep = "\t", index = False)

