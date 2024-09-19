## This utility script roots a tree using a supplied outgroup taxon, clips the outgroup from the tree and saves it at the supplied location.
##
## USAGE
## Rscript remove_outgroup.R @tree @output @outgroup_label
##
## PARAMETERS
## tree             a phylogeny in Newick format
## output           path to the directory in which the rooted tree will be saved
## outgroup_label   label of the outgroup taxon in the tree

library(ape)

args = commandArgs(trailingOnly=TRUE)
tree.file = args[1]
out.dir = args[2]
outgroup = args[3]

tree.0 =  read.tree(tree.file)
tree.m = root(tree.0, outgroup = outgroup)
tree.m = drop.tip(tree.m, tip = outgroup)

write.tree(tree.m, file = paste(out.dir, "merge.contree", sep = "/"))
