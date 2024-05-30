#' Preprocess SingleCellExperiment

library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(yaml)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Preprocess SingleCellExperiment")
parser$add_argument('--sce_rds_input', metavar='DIR', type='character',
                    help="Path to input rds")
parser$add_argument('--metadata', metavar='FILE', type='character',
                    help="Path to metadata yaml file")
parser$add_argument('--sample', type='character',
                    help="Sample ID")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path.")
args <- parser$parse_args()

sce <- readRDS(args$sce_rds_input)

# Read metadata
metadata <- read_yaml(args$metadata)
metadata_df <- dplyr::bind_rows(lapply(names(metadata), function(i) {
  data.frame(metadata[[i]])
}))
print(metadata_df)
print(args$sample)
metadata_row <- metadata_df %>%
  dplyr::filter(id == args$sample)
print("or here")
# Add metadata
print(metadata_row)
colData(sce) <- colData(sce) %>%
  cbind(metadata_row)
print("and here")
mito_genes <- str_detect(rowData(sce)$Symbol, "^MT\\-")

ribo_genes <- str_detect(rowData(sce)$Symbol, "^RP(L|S)")
print("here")
# Calculate basic QC stats 
sce <- calculateQCMetrics(sce, exprs_values = "counts", feature_controls =
                            list(mitochondrial=mito_genes, ribosomal=ribo_genes))

# Write outputs
saveRDS(sce, file = args$outfname)

cat("Completed.\n")

