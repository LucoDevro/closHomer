#!/usr/bin/env python3

import pandas as pd
import sys

def basename(path):
    return path.split('/')[-1].split('.')[0]

sec_clust_df_file = sys.argv[1]

sec_clust_df = pd.read_table(sec_clust_df_file, converters={0: basename, 1: basename})

cluster_sizes = sec_clust_df.groupby("nearest_representative_genome")['genome'].count()
average_ani = sec_clust_df.groupby("nearest_representative_genome")['average_nucleotide_identity'].mean()
average_af = sec_clust_df.groupby("nearest_representative_genome")['alignment_fraction'].mean()

cluster_metrics = pd.concat([cluster_sizes, average_ani, average_af], axis = 1)
cluster_metrics.to_csv('cluster_metrics', sep = "\t")