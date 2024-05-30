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
parser$add_argument('--mito', type='double',
                    help="Mitochondrial threshold", default = 20)
parser$add_argument('--ribo', type='double',
                    help="Ribosomal threshold", default = 60)
parser$add_argument('--nmads', type='double',
                    help="Number of MADs to filter at.", default = 3)
parser$add_argument('--features', type='double',
                    help="Minimum number of features", default = 1000)
parser$add_argument('--min_malat', type='double',
                    help="Minimum MALAT1 expression", default = 0)
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for preprocessed SCE.")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

# Filter SCE
sce_filtered <- filter_cells(sce, nmads = args$nmads, type = "lower", 
                             log = TRUE, max_mito = args$mito, max_ribo = args$ribo)

# save.image('/cellassign/fitness-scrna/results/dumps/filter-dump.Rdata')

sce_filtered <- sce_filtered[, sce_filtered$total_features_by_counts >= args$features]

# MA: for the already normalized files we can remove the row below, as we already have rownames
# for the old unnormalized files we need this, but for the new ones we don't
rownames(sce_filtered) <- scater::uniquifyFeatureNames(rowData(sce_filtered)$ID, rowData(sce_filtered)$Symbol)

if ("logcounts" %in% names(assays(sce_filtered))) {
  sce_filtered <- sce_filtered[, as.vector(logcounts(sce_filtered)['MALAT1',]) > args$min_malat]
} else {
  sce_filtered <- sce_filtered[, as.vector(log2(counts(sce_filtered)['MALAT1',])) > args$min_malat]
}

# Write outputs
saveRDS(sce_filtered, file = args$outfname)

cat("Completed.\n")

