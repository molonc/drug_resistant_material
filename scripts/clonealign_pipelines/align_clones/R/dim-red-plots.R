
suppressPackageStartupMessages({
  library(scater)
  library(tidyverse)
  library(clonealign)
})

sce <- readRDS("/cellassign/fitness-scrna/results/v5/outputs/preprocess/sce_annotated/SA906_p15b.rds")
ca <- readRDS("/cellassign/fitness-scrna/results/v5/outputs/align_clones/clonealign_fit/SA906_p15b.rds")

ca <- clonealign:::recompute_clone_assignment(ca, 0.5)

colnames(sce) <- sce$Barcode

barcodes <- ca$clone_fit$Barcode

sce <- sce[, barcodes]
sce$clone <- ca$clone

plotTSNE(sce, colour_by = "clone")
