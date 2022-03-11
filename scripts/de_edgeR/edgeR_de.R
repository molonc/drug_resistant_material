suppressPackageStartupMessages({
  library(knitr)
  library(tidyverse)
  # library(tensorflow)
  library(SingleCellExperiment)
  library(scater)
  library(data.table)
  # library(pheatmap)
  # library(ggrepel)
  library(grid)
  library(scran)
  library(yaml)
  library(scales)
  # library(kableExtra)
  # library(dendextend)
  library(edgeR)
  library(org.Hs.eg.db)
  library(Seurat)
  # library(pheatmap)
  library(clonealign)
  # library(ggrepel)
  library(annotables)
  library(scales)
  library(RColorBrewer)
  # library(fgsea)
  
  # library(scrna.utils)
  # library(scrna.sceutils)
  # library(cellassign)
  # library(cellassign.utils)
  
  # library(svglite)
})


source("/home/htran/Projects/farhia_project/rscript/de_edgeR/de-utils.R")

load_data_edgeR <- function(sce_fn, clonelabel1, clonelabel2,
                            save_dir, datatag){
  
  output_dir <- paste0(save_dir, datatag)
  # output_fig_dir <- paste0("../results/", params$series, "-", params$version, "/de")
  dir.create(output_dir, recursive=TRUE)
  # 
  # ca_dir <- paste0("../../results/", params$series, "-", params$version, "/outputs/align_clones/clonealign_fit_csv/")  
  # 
  # clone_dir <- paste0("../../results/", params$series, "-", params$version, "/outputs/parse_cnv/gene_clone_cn/")
  
  # sce_dir <- paste0("../../results/", params$series, "-", params$version, "/outputs/preprocess/sce_twice_scran/")
  # title <- params$title
  # input_sce1 <- paste0(sce_dir, params$sample1, ".rdata")
  # input_ca1 <- paste0(ca_dir, params$sample1, ".csv")
  # clone_cn1 <- paste0(clone_dir, params$sample1, ".tsv")
  # clone1 <- params$clone1
  # clonelabel1 <- params$clonelabel1
  
  
  # input_sce2 <- paste0(sce_dir, params$sample2, ".rdata")
  # input_ca2 <- paste0(ca_dir, params$sample2, ".csv")
  # clone_cn2 <- paste0(clone_dir, params$sample2, ".tsv")
  # clone2 <- params$clone2
  # clonelabel2 <- params$clonelabel2
  
  # Sometimes clone1 or clone2 can be for example B_R, then look for all the clones in dlp
  # clone1list <- unlist(strsplit(clone1, "_"))
  # clone2list <- unlist(strsplit(clone2, "_"))
  
  # sample <- paste0(params$cno,"-", params$sample1,"-", clone1,"-", params$sample2, "-", clone2)
  # 
  # print(output_fig_dir)
  # output_file <- glue("{output_fig_dir}/{sample}_summary.html")
  
  
  
  # sce1 <- readRDS(input_sce1)
  # sce1 <- fixsce(sce1)
  # # this correction should be done in identify-murine.R
  # sce1 <- sce1[which(!is.na(rowData(sce1)$ID)),]
  # df_cnv1 <- read_tsv(clone_cn1)
  # df_cnv1 <- df_cnv1[df_cnv1$cluster %in% clone1list,]
  # df_cnv1$cluster <- clonelabel1
  # colnames(sce1) <- sce1$Barcode
  #ca1 <- readRDS(input_ca1)
  #ca1 <- clonealign:::recompute_clone_assignment(ca1, 0.5)
  #sce1$clone <- ca1$clone
  #sce_qc1 <- sce1[, sce1$clone %in% clone1]
  #mybarcodes1 <- ca1$clone_fit[ca1$clone_fit$clone==clone1,"Barcode"]
  # ca1 <- read.csv(input_ca1,sep=",")
  # mybarcodes1 <- ca1[ca1$clone %in% clone1list,"Barcode"]
  # sce_qc1 <- sce1[,colData(sce1)$Barcode %in% mybarcodes1]
  ### MA 23 Apr 2020: I used to have 1_, but 2_ will actually put it first
  sce_qc1$clone <- paste0("2_", clone1)
  # remove mitocondrial genes and other genes that are not of interest
  rowData(sce_qc1)$Symbol <- rownames(counts(sce_qc1))
  #sce_qc1 <- sce_qc1[!grepl("^RP[L|S]|^MT-|^FOS|^JUN|^HSP", rowData(sce_qc1)$Symbol)]
  
  sce2 <- readRDS(input_sce2)
  sce2 <- fixsce(sce2)
  sce2 <- sce2[which(!is.na(rowData(sce2)$ID)),]
  df_cnv2 <- read_tsv(clone_cn2)
  df_cnv2 <- df_cnv2[df_cnv2$cluster %in% clone2list,]
  df_cnv2$cluster <- clonelabel2
  colnames(sce2) <- sce2$Barcode
  #ca2 <- readRDS(input_ca2)
  #ca2 <- clonealign:::recompute_clone_assignment(ca2, 0.5)
  #mybarcodes2 <- ca2$clone_fit[ca2$clone_fit$clone==clone2,"Barcode"]
  ca2 <- read.csv(input_ca2,sep=",")
  mybarcodes2 <- ca2[ca2$clone %in% clone2list,"Barcode"]
  sce_qc2 <- sce2[,colData(sce2)$Barcode %in% mybarcodes2]
  sce_qc2$clone <- paste0("1_", clone2)
  # remove mitocondrial genes and other genes that are not of interest
  rowData(sce_qc2)$Symbol <- rownames(counts(sce_qc2))
  #sce_qc2 <- sce_qc2[!grepl("^RP[L|S]|^MT-|^FOS|^JUN|^HSP", rowData(sce_qc2)$Symbol)]
  
  if (sum(rowData(sce1)$ID!=rowData(sce2)$ID) > 0) {
    stop("Ensemble gene names are different between the two sce samples!")
  }
  
  #ensgenes <- rowData(sce1)$ID
  
  # combine the two sces
  saveRowData <- rowData(sce_qc1)
  rowData(sce_qc1) <- rowData(sce_qc2) <- NULL
  sce_qc <- cbind(sce_qc1, sce_qc2)
  # Note: keeping the data for the first one, so mean_counts, log10_mean_counts will not be correct for the bind
  rowData(sce_qc) <- saveRowData
  df_cnv <- rbind(df_cnv1, df_cnv2)
  sce_qc <- sce_qc[, sce_qc$total_features_by_counts > 1500]
}

edgeR_de_analysis <- function(){
  sce_de <- sce_qc
  
  rs <- rowSums(as.matrix(counts(sce_qc)))
  qplot(rs, log='x') + geom_vline(xintercept = 100)
  
  sce_de <- sce_de[rowSums(as.matrix(counts(sce_de))) > 100, ]
  mycounts <- as.matrix(counts(sce_de))    # zeros are good
  dge <- DGEList(counts=mycounts, group=sce_de$clone)
  
  # which vs which ???
  design <- model.matrix(~ clone, data = colData(sce_de))
  
  # This describes the edgeR user manual
  # http://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
  
  dge <- estimateDisp(dge, design = design)
  fit <- glmQLFit(dge, design = design)
  qlf <- glmQLFTest(fit)
  tt <- topTags(qlf, n = Inf)
  
  tt <- as.data.frame(tt) %>% 
    rownames_to_column("gene_symbol")
  tt$ensembl_gene_id <- rowData(sce_de)$ID[match(tt$gene_symbol,rowData(sce_de)$Symbol)]
  
  # Saving the logFC file
  ttfile <- glue("{output_fig_dir}/{sample}_logfc_results.csv")
  write.table(tt,file=ttfile,sep=",", row.names=FALSE, quote=FALSE)
  # tt_df_fdr <- tt_df[tt_df$FDR<0.01,]
}

