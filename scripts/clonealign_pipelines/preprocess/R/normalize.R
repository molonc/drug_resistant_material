#' Preprocess SingleCellExperiment

library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Preprocess SingleCellExperiment")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--umap_neighbors', type='integer',
                    help="Number of nearest neighbours for UMAP", default = 15)
parser$add_argument('--umap_min_dist', type = 'double',
                    help="Minimum distance for UMAP", default = 0.1)
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for preprocessed SCE.")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

# Remove doublets
if ("doublet" %in% colnames(colData(sce))) {
  sce <- sce %>%
    scater::filter(doublet == "Singlet")
}

# Compute size factors
qclust <- quickCluster(sce, min.size = 100)
sce <- computeSumFactors(sce, clusters = qclust)

sce$size_factor <- sizeFactors(sce)

# Normalize
# TODO: Implement other methods for normalization (e.g. see Rafa Iziarry paper)
# MA: Not sure if I should normalize here, I am removing it for now
# 
sce_normalized <- sce

# Dimensionality reduction
if ("logcounts" %in% names(assays(sce_normalized))) {
  sce_normalized <- runPCA(sce_normalized, ntop = 1000, ncomponents = 50, exprs_values = "logcounts")
  sce_normalized <- runTSNE(sce_normalized, use_dimred = "PCA", n_dimred = 50, ncomponents = 2)
  sce_normalized <- runUMAP(sce_normalized, use_dimred = "PCA", n_dimred = 50, ncomponents = 2, n_neighbors = args$umap_neighbors, min_dist = args$umap_min_dist)
}
# Write outputs
saveRDS(sce_normalized, file = args$outfname)

cat("Completed.\n")



