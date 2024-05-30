#' Cell cycle prediction with cyclone on SingleCellExperiment

library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(org.Hs.eg.db)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Run cyclone on SingleCellExperiment")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--ncpus', type='double',
                    help="Number of cores to use.")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for cell cycle assignments.")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))

## MA: There is a 1 to many correspondence 
if(is.null(rowData(sce)$ID)) {
  rowData(sce)$ID <- rownames(sce)
}

# Run cyclone
if (!str_detect(rowData(sce)$ID, "^ENSG")) {
  gene_ids <- mapIds(org.Hs.eg.db, rowData(sce)$ID, 'ENSEMBL', 'SYMBOL')
  sce <- sce[!is.na(gene_ids),]
  gene_ids <- unname(gene_ids[!is.na(gene_ids)][rowData(sce)$ID])
} else {
  gene_ids <- rowData(sce)$ID
}

assignments <- cyclone(sce, hs.pairs, gene.names=gene_ids, 
                       min.iter = 10, verbose = TRUE, BPPARAM = MulticoreParam(args$ncpus))

# Write outputs
saveRDS(assignments, file = args$outfname)

cat("Completed.\n")


