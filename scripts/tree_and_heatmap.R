# Installing packages.
BiocManager::install("pheatmap")
BiocManager::install("data.table")

# Importing packages.
library(tidyr)
library(ape)
library(ggplot2)
library(pheatmap)
library(glue)
library(data.table)
library(RColorBrewer)

# I/O paths.
# cutoff <- snakemake@params[[1]]
# mashtable <- snakemake@input[[1]]
# metadata <- snakemake@input[[2]]
# output <- snakemake@output[[1]]

### Testing paths. ###
cutoff <- 0.01
mashtable <-  '/Users/januszkoszucki/Work/BPdiversity/input/distances.tab'
metadata <- '/Users/januszkoszucki/Work/BPdiversity/output/metadata.csv'
phyloheatmap.out <- glue('/Users/januszkoszucki/Work/BPdiversity/output/test-phyloheatmap2-{cutoff}.pdf')

# Prepare vector colors for K locus as category (max 74 categories!).
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
colors_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

##############################################
############### MAIN #########################
##############################################

# Loading distance (mash) table.
dist.df <- read.table(mashtable, header=FALSE, sep="")
dist.df <- dist.df[1:3]
colnames(dist.df) <- c("query", "subject", "mash")
dm <- spread(dist.df, 'subject', 'mash')

# Loading metadata (ST, KL, phage variants, in future: phage RBPs).
metadata.df <- read.csv(metadata, sep=',', row.names=1)

# Clearing distance (mash) table.
rownames(dm) <- dm[,'query']
dm <- dm[,2:length(dm)]
labels <- rownames(dm)

# Get genomes representative genomes (by indicated cutoff).
dm.dist <- as.dist(dm)
lineages <- hclust(dm.dist)
lineages.clusters <- cutree(lineages, h=cutoff)
representative.values <- unique(lineages.clusters)
representative.idx <- match(representative.values, lineages.clusters)
representative.genomes <- labels[representative.idx]

# Adjust metadata matrix to have only rows and columns that has been chosen as representative genomes from hclust algorithm.
dm <- dm[representative.genomes,representative.genomes]
metadata.df <- metadata.df[representative.genomes,]
heatmap.variants <- metadata.df[,4:length(metadata.df)]
# Remove "PV" string from column names of phage variants.
colnames(heatmap.variants) <- sub('PV', '', colnames(heatmap.variants))
# Remove phage variants that are not present in selected representative genomes.
heatmap.variants <- heatmap.variants[,colSums(heatmap.variants[-1,])>1]


# Extract to separate data frame information about K locus types to annotate rows.
klocus.df <- data.frame(metadata.df[representative.genomes, 'K_locus'])
rownames(klocus.df) <- rownames(metadata.df)
colnames(klocus.df) <- c('K locus')


# Plot heatmap.
plot.title <- glue('phageVariants vs klebsiella (cutoff: {cutoff})')
ncategories <- length(unique(klocus.df[['K locus']]))
plot <- pheatmap(heatmap.variants, annotation_row = klocus.df, clustering_distance_rows=as.dist(dm),border_color='black', fontsize_row = 10, fontsize_col = 14, main=plot.title, annotation_colors=colors_vector[1:ncategories])

# Save heatmap.
save_pheatmap_pdf <- function(x, filename, width=16, height=16) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(plot, output)


####################################
############ Next step #############
####################################

### Plot all capsules vs all phages from a cluster accordingly to phylogenetic tree ###



# # Reshape data frame to pivot (binary) table.
# ncategories <- length(unique(klocus.df[['K locus']]))
# klocus.df[['values']] <- 1
# klocus.df[['ID']] <- rownames(klocus.df)
# klocus.pivot <- spread(klocus.df, 'K locus', 'values')
# klocus.pivot[is.na(klocus.pivot)] <- 0
# rownames(klocus.pivot) <- rownames(klocus.df)
# klocus.pivot <- klocus.pivot[,2:ncategories+1]