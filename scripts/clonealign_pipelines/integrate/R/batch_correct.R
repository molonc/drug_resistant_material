#' Integrate samples into groups

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

parser <- ArgumentParser(description = "Integrate SCEs")
parser$add_argument('--sce', metavar='FILE', type='character', nargs = '+',
                    help="Path(s) to SCEs")
parser$add_argument('--batch_variable', type='character',
                    help="Batch variable")
parser$add_argument('--method', type='character',
                    help="Integration method", default = 'scanorama')
parser$add_argument('--rowdata_cols_remove', type='character', nargs = '+',
                    help="Columns to remove in rowdata", default = c("mean_counts", "log10_mean_counts", "n_cells_by_counts", "pct_dropout_by_counts", "total_counts", "log10_total_counts"))
parser$add_argument('--conda_env', type = 'character',
                    help="Name of conda environment with tensorflow", default = "r-tensorflow")
parser$add_argument('--conda_path', type = 'character',
                    help="Conda path", default = "auto")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for preprocessed SCE.")
args <- parser$parse_args()

# Reset snakemake default
Sys.setenv(PYTHONPATH='')

reticulate::use_condaenv(args$conda_env, 
                         conda = args$conda_path)

sce_paths <- unlist(args$sce)
sces <- lapply(sce_paths, function(x) {
  sce <- readRDS(x)
  keep_cols <- which(!colnames(rowData(sce)) %in% unlist(args$rowdata_cols_remove))
  rowData(sce) <- rowData(sce)[,keep_cols]
  
  rowData(sce)$start <- NULL
  rowData(sce)$end <- NULL
  rowData(sce)$strand <- NULL
  
  return(sce)
})

# save.image(paste0("/cellassign/fitness-scrna/results/dumps/_integrate-dump-", 1, ".Rdata"))

message("Binding SCEs ...")
sce_merged <- do.call('cbind', sces)

# Recompute rowdata QC
sce_merged <- calculateQCMetrics(sce_merged)

# Method for batch correction
batch_correct <- function(sce, batch_col, method = "scanorama") {
  ## Feature selection
  batches <- unique(colData(sce)[,batch_col])
  
  fit <- trendVar(sce, parametric=TRUE, use.spikes = FALSE)
  decomp <- decomposeVar(sce, fit)
  decomp$Symbol <- rowData(sce)$Symbol
  
  decomp_chosen <- decomp %>% subset(bio > 0)
  chosen <- rownames(decomp_chosen)
  
  chosen_features <- rep(list(chosen), length(batches))
  
  ## Batch correction
  if (method == "scanorama") {
    norm_matrices <- lapply(batches, function(bat) {
      sce_sub <- sce[,colData(sce)[,batch_col] == bat]
      norm_data <- t(as.matrix(logcounts(sce_sub)))
      
      # Discovered on 2019/05/06: scanorama DOES NOT properly subset
      norm_data <- norm_data[,chosen]
      return(norm_data)
    })
    
    indexes <- lapply(batches, function(bat) {
      which(colData(sce)[,batch_col] == bat)
    })
    
    scanorama <- reticulate::import('scanorama')
    
    message("Running scanorama ...")
    # Integration and batch correction
    integrated_corrected_data <- scanorama$correct(norm_matrices,
                                                   chosen_features,
                                                   return_dimred = TRUE,
                                                   return_dense = TRUE)
    
    scanorama_mat <- matrix(NA, nrow = sum(sapply(indexes, length)), ncol = 100)
    scanorama_mat_expr <- matrix(NA, nrow = sum(sapply(indexes, length)), ncol = ncol(integrated_corrected_data[[2]][[1]]))
    for (i in seq_along(indexes)) {
      idx <- indexes[[i]]
      scanorama_mat[idx,] <- integrated_corrected_data[[1]][[i]]
      scanorama_mat_expr[idx,] <- integrated_corrected_data[[2]][[i]]
    }
    reducedDim(sce, "scanorama_int") <- scanorama_mat
    
    sce_sel <- sce
    sce_sel <- sce_sel[as.character(integrated_corrected_data[[3]]),]
    assay(sce_sel, "scanorama") <- t(scanorama_mat_expr)
    pca_res <- pca2(sce_sel, ntop = 1000, ncomponents = 50, exprs_values = "scanorama")
    reducedDim(sce_sel, "PCA2") <- pca_res$x
    
    sce_sel <- runTSNE(sce_sel, use_dimred = "PCA2", n_dimred = 50)
    sce_sel <- runUMAP(sce_sel, use_dimred = "PCA2", n_dimred = 50)
    
    reducedDim(sce, "scanorama_PCA") <- reducedDim(sce_sel, "PCA2")
    reducedDim(sce, "scanorama_TSNE") <- reducedDim(sce_sel, "TSNE")
    reducedDim(sce, "scanorama_UMAP") <- reducedDim(sce_sel, "UMAP")
    
    colnames(scanorama_mat_expr) <- integrated_corrected_data[[3]]
    reducedDim(sce, "scanorama_exprs") <- scanorama_mat_expr
    
    sce@metadata$batchcor_pca_res <- pca_res
    sce@metadata$batchcor_genes <- chosen
  } else if (method == "harmony") {
    library(harmony)
    
    harmony_embeddings <- HarmonyMatrix(logcounts(sce)[chosen,], colData(sce) %>% as.data.frame, batch_col, do_pca = TRUE)
    colnames(harmony_embeddings) <- paste0("PC", 1:ncol(harmony_embeddings))
    
    sce2 <- sce
    reducedDim(sce2, "harmony_PCA") <- harmony_embeddings
    
    sce2 <- runTSNE(sce2, use_dimred = "harmony_PCA", n_dimred = ncol(harmony_embeddings))
    sce2 <- runUMAP(sce2, use_dimred = "harmony_PCA", n_dimred = ncol(harmony_embeddings))
    
    reducedDim(sce, "harmony_PCA") <- reducedDim(sce2, "harmony_PCA")
    reducedDim(sce, "harmony_TSNE") <- reducedDim(sce2, "TSNE")
    reducedDim(sce, "harmony_UMAP") <- reducedDim(sce2, "UMAP")
    sce@metadata$batchcor_genes <- chosen
  } else {
    stop("Other methods not implemented.")
  }
  
  return(sce)
}

message("Running batch correction ...")
sce_bc <- batch_correct(sce_merged, batch_col = args$batch_variable)

# Write outputs
saveRDS(sce_bc, file = args$outfname)

cat("Completed.\n")

