## This utility script extracts a subset of taxa listed in a given index file from a provided phylogeny, and saves at the new subtree at the supplied path.
##
## USAGE
## Rscript extract_tree.R @tree @index @out
##
## PARAMS
## tree     a phylogeny in Newick format
## index    a text file containing all taxa to be extracted, separated by newlines
## out      path to where the new subtree should be saved

args = commandArgs(trailingOnly=TRUE)
tree.file = args[1]
index = args[2]
out = args[3]

library(ape)

tree.m = read.tree(tree.file)
tips = as.character(read.table(index)[,1])
tree = keep.tip(tree.m, tip = tips)
write.tree(tree, file = out)
