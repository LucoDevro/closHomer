## This utility script plots an MDS plot of the Mash pairwise distances between all samples and adds the species clusters clusters as an additional metadata colour scale.
##
## USAGE
## Rscript plot_mash_grouped.R @mds_coords @metadata @out_dir
## mds_coords   tsv file with the 2D coordinates of the MDS representation of all samples as produced by panaroo-qc
## metadata     tsv file with the metadata of the samples, should at least contain a column "Genome_accession" with the RefSeq Assembly identifiers of the samples, and a "Taxonomic_cluster" column with the species cluster assignment
## out_dir      output directory for the MDS plot

args = commandArgs(trailingOnly=TRUE)
mds_coords.file = args[1]
metadata.file = args[2]
output.folder = args[3]

coords = read.table(mds_coords.file, sep='\t', header = T)
colnames(coords)[1] = 'Genome_accession'
metadata = read.table(metadata.file, sep='\t', header = T)

svg(paste(output.folder, 'MDS_mash_plot_clusters.svg', sep = "/"))
coords.with.cluster = subset(merge(coords, metadata, by = "Genome_accession"), select = c('Genome_accession','coordx','coordy','Taxonomic_cluster'))
plot(coords.with.cluster$coordx, coords.with.cluster$coordy, pch = 19, col = factor(coords.with.cluster$Taxonomic_cluster), xlab='MDS1', ylab='MDS2')
legend('topleft', legend = levels(factor(coords.with.cluster$Taxonomic_cluster)), pch=19, col = factor(levels(factor(coords.with.cluster$Taxonomic_cluster))))
dev.off()
