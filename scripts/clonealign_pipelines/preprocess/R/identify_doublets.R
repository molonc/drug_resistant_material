#' Identify doublets in SingleCellExperiment

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


parser <- ArgumentParser(description = "Remove doublets from SingleCellExperiment")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--method', type='character',
                    help="Doublet method", default = 'scrublet')
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

sce_path <- args$sce
sce <- readRDS(sce_path)

sce <- infer_doublets(sce, 
                      method = args$method,
                      conda_binary = "/home/rstudio/miniconda/bin/conda")

# Write outputs
saveRDS(sce, file = args$outfname)

cat("Completed.\n")



### function from cellassign.utils
#' #' Infer doublets
#' #'
#' #' @param sce SingleCellExperiment
#' #' @param group_col Group column (for DoubletFinder)
#' #' @param method Doublet method
#' #' @param python_path Python path to use
#' #' @param conda_env Conda environment to use
#' #' @param conda_binary Conda binary path to use
#' #' @param min_total_genes Minimum total genes to use as markers (for DoubletDecon)
#' #' @export
#' infer_doublets <- function(sce,
#'                            group_col = NULL,
#'                            method = "doubletfinder",
#'                            python_path = NULL,
#'                            conda_env = "r-tensorflow",
#'                            conda_binary = "auto",
#'                            min_total_genes = 50) {
#'   if (!is.null(group_col)) {
#'     if (!group_col %in% colnames(colData(sce))) {
#'       stop("Group column specified does not exist!")
#'     }
#'   }
#'   if (method %in% c("scrublet", "doubletdetection")) {
#'     if (!is.null(python_path)) {
#'       reticulate::use_python(python_path, required = FALSE)
#'     } else {
#'       if (!is.null(conda_env)) {
#'         reticulate::use_condaenv(conda_env, conda = conda_binary,
#'                                  required = FALSE)
#'       } else {
#'         stop("Environment could not be started.")
#'       }
#'     }
#'   }
#'   if (method == "doubletfinder") {
#'     library(Seurat)
#'     library(DoubletFinder)
#'     library(fields)
#'     library(modes)
#'     sce2 <- sce
#'     colnames(sce2) <- paste0(sce2$Barcode, "_", as.numeric(factor(sce2$Sample)))
#'     seurat_obj <- Convert(from = sce2, to = "seurat",
#'                           min.genes = -1)
#'     seurat_obj@dr <- list()
#'     seurat_obj <- NormalizeData(seurat_obj)
#'     seurat_obj <- ScaleData(seurat_obj, vars.to.regress = "total_counts")
#'     seurat_obj <- FindVariableGenes(seurat_obj, x.low.cutoff = 0.0125, y.cutoff = 0.25, do.plot=FALSE)
#'     seurat_obj <- RunPCA(seurat_obj, pc.genes = seurat_obj@var.genes, pcs.print = 0)
#'     seurat_obj <- RunTSNE(seurat_obj, dims.use = 1:10, verbose=TRUE)
#'     seurat_obj <- FindClusters(object = seurat_obj, reduction.type = "pca", dims.use = 1:10,
#'                                resolution = 0.6, print.output = 0, save.SNN = TRUE)
#'     sweep.res.list <- DoubletFinder:::paramSweep(seurat_obj)
#'     sweep.stats <- DoubletFinder:::summarizeSweep(sweep.res.list, GT = FALSE)
#'     bcmvn <- DoubletFinder:::find.pK(sweep.stats)
#'     pk_summary <- sweep.stats %>%
#'       dplyr::group_by(pK) %>%
#'       dplyr::summarise(BCsd=sd(BCreal), BCreal=mean(BCreal), BCmetric=BCreal/BCsd^2) %>%
#'       dplyr::ungroup()
#'     selected_pk <- as.numeric(as.character(pk_summary$pK[which.max(pk_summary$BCmetric)]))
#'     if (is.null(group_col)) {
#'       seurat_obj <- FindClusters(object = seurat_obj, reduction.type = "pca", dims.use = 1:20,
#'                                  resolution = 0.6, print.output = 0, save.SNN = TRUE)
#'       seurat_obj@meta.data$seurat_cluster <- seurat_obj@ident
#'       group_col <- "seurat_cluster"
#'     } else {
#'       seurat_obj@ident <- factor(seurat_obj@meta.data[,group_col])
#'     }
#'     homotypic.prop <- DoubletFinder:::modelHomotypic(seurat_obj@meta.data[,group_col])
#'     nExp_poi <- round(0.075*length(seurat_obj@cell.names))
#'     nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#'     seurat_obj <- DoubletFinder:::doubletFinder(seurat_obj, pN = 0.25, pK = selected_pk, nExp = nExp_poi, reuse.pANN = FALSE)
#'     seurat_obj <- DoubletFinder:::doubletFinder(seurat_obj,
#'                                                 pN = 0.25,
#'                                                 pK = selected_pk,
#'                                                 nExp = nExp_poi.adj,
#'                                                 reuse.pANN = colnames(seurat_obj@meta.data)[str_detect(colnames(seurat_obj@meta.data),
#'                                                                                                        paste0("pANN_", 0.25, "_",
#'                                                                                                               selected_pk))][1])
#'     DF_classification_cols <- colnames(seurat_obj@meta.data)[str_detect(colnames(seurat_obj@meta.data),
#'                                                                         "^DF\\.classifications")]
#'     pANN_cols <- colnames(seurat_obj@meta.data)[str_detect(colnames(seurat_obj@meta.data),
#'                                                            "^pANN")]
#'     seurat_obj@meta.data$DF_doublet <- ifelse((seurat_obj@meta.data[,DF_classification_cols[1]] == "Doublet" & seurat_obj@meta.data[,DF_classification_cols[2]] == "Doublet"),
#'                                               "Doublet_High",
#'                                               ifelse((seurat_obj@meta.data[,DF_classification_cols[1]] == "Singlet" & seurat_obj@meta.data[,DF_classification_cols[2]] == "Singlet"),
#'                                                      "Singlet",
#'                                                      "Doublet_Low"))
#'     sce$DF_doublet <- seurat_obj@meta.data$DF_doublet
#'     sce$doublet <- sce$DF_doublet
#'   } else if (method == "doubletdetection") {
#'     dd <- reticulate::import(module = 'doubletdetection')
#'     raw_counts <- as.matrix(Matrix::t(counts(sce)))
#'     clf <- dd$BoostClassifier()
#'     labels <- clf$fit(raw_counts)$predict()
#'     sce$DD_doublet <- labels %>%
#'       plyr::mapvalues(from = c(0, 1),
#'                       to = c("Singlet", "Doublet"))
#'     sce$doublet <- sce$DD_doublet
#'   } else if (method == "scrublet") {
#'     scr <- reticulate::import(module = 'scrublet')
#'     raw_counts <- as.matrix(Matrix::t(counts(sce)))
#'     scrub <- scr$Scrublet(raw_counts, expected_doublet_rate = 0.04)
#'     scrub_result <- scrub$scrub_doublets()
#'     sce$scrublet_scores <- scrub_result[[1]]
#'     sce$scrublet_class <- scrub_result[[2]] %>%
#'       plyr::mapvalues(c(FALSE, TRUE),
#'                       c("Singlet", "Doublet"))
#'     sce$doublet <- sce$scrublet_class
#'   } else if (method == "doubletdecon") {
#'     library(DoubletDecon)
#'     library(Seurat)
#'     sce2 <- sce
#'     colnames(sce2) <- paste0(sce2$Barcode, "_", as.numeric(factor(sce2$Sample)))
#'     seurat_obj <- Convert(from = sce2, to = "seurat",
#'                           min.genes = -1)
#'     seurat_obj@dr <- list()
#'     seurat_obj <- NormalizeData(seurat_obj)
#'     seurat_obj <- ScaleData(seurat_obj, vars.to.regress = "total_counts")
#'     seurat_obj <- FindVariableGenes(seurat_obj, x.low.cutoff = 0.0125, y.cutoff = 0.25, do.plot=FALSE)
#'     seurat_obj <- RunPCA(seurat_obj, pc.genes = seurat_obj@var.genes, pcs.print = 0)
#'     if (is.null(group_col)) {
#'       seurat_obj <- FindClusters(object = seurat_obj, reduction.type = "pca", dims.use = 1:20,
#'                                  resolution = 0.6, print.output = 0, save.SNN = TRUE)
#'       seurat_obj@meta.data$seurat_cluster <- seurat_obj@ident
#'     } else {
#'       groups <- factor(as.numeric(factor(seurat_obj@meta.data[,group_col])) - 1)
#'       names(groups) <- seurat_obj@cell.names
#'       seurat_obj@ident <- groups
#'     }
#'     seurat_markers <- FindAllMarkers(object = seurat_obj, only.pos = TRUE, min.pct = 0.25,
#'                                      thresh.use = 0.25)
#'     num_genes_per_cluster <- max(ceiling(min_total_genes/length(unique(seurat_markers$cluster))), 5)
#'     seurat_top_markers <- seurat_markers %>%
#'       dplyr::group_by(cluster) %>%
#'       dplyr::top_n(num_genes_per_cluster, avg_logFC)
#'     temp_out_dir <- tempdir()
#'     expression_path <- file.path(temp_out_dir, "expression_mat.txt")
#'     marker_gene_path <- file.path(temp_out_dir, "top_genes.txt")
#'     cluster_path <- file.path(temp_out_dir, "cluster.txt")
#'     raw_counts <- as.matrix(counts(sce2))
#'     write.table(raw_counts, file = expression_path, row.names = TRUE, col.names = NA, sep = "\t")
#'     write.table(seurat_top_markers, file = marker_gene_path, row.names = TRUE, col.names = NA, sep = "\t")
#'     write.table(seurat_obj@ident, file = cluster_path, row.names = TRUE, col.names = NA, sep = "\t")
#'     newfiles <- Seurat_Pre_Process(expression_path, marker_gene_path, cluster_path)
#'     newfiles$newGroupsFile <- as.data.frame(newfiles$newGroupsFile)
#'     results <- Main_Doublet_Decon(rawDataFile=newfiles$newExpressionFile,
#'                                   groupsFile=newfiles$newGroupsFile,
#'                                   filename="test123",
#'                                   location=temp_out_dir,
#'                                   fullDataFile=NULL,
#'                                   removeCC=FALSE,
#'                                   species="hsa",
#'                                   rhop=1,
#'                                   write=TRUE,
#'                                   PMF=TRUE,
#'                                   useFull=FALSE,
#'                                   heatmap=FALSE,
#'                                   centroids=FALSE,
#'                                   num_doubs=100,
#'                                   downsample="none",
#'                                   sample_num=NULL,
#'                                   only50=FALSE,
#'                                   min_uniq=4)
#'   } else {
#'     stop("Unrecognized method.")
#'   }
#'   return(sce)
#' }



