suppressPackageStartupMessages({
  require("optparse")
  require("scater")
  require("argparse")
  require("SingleCellExperiment")
  require("scMerge")
  require("dplyr")
  require("ggplot2")
  require("stringr")
  
})




plot_variation_function <- function(meta_data, xstring="treatment_status", ystring="gene_exp", plottype="batches_info", plottitle="Raw data",
                                    xlabel='', ylabel="Stably express genes") {
  
  p <- ggplot(meta_data, aes_string(x=xstring, y=ystring, fill=plottype)) +
    stat_boxplot(aes_string(xstring, ystring), geom='errorbar', linetype=1, width=0.5)+ #whiskers
    geom_boxplot(aes_string(xstring, ystring),outlier.shape=1) +
    stat_summary(fun.y=mean, geom="point", size=2) +
    stat_summary(fun.data = mean_se, geom = "errorbar")
  
  p <- p + labs(x=xlabel,y=ylabel,title=plottitle)
  p <- p + theme(legend.title = element_text(color="black", size=11, hjust = 0.5),
                 plot.title = element_text(color="black", size=14, hjust = 0.5),
                 # legend.position = "none",
                 # axis.line = element_blank(),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 axis.text.x = element_text(angle = 90, hjust = 1,size=9),
                 #axis.title = element_blank(),
                 axis.ticks = element_blank()
  )
  # p
  return(p)
  
}

get_plt_SEG <- function(sce, seglist, plottitle='Batch effect estimation', yl=c(0,1.6)){
  sce_SEG <- sce[seglist,]
  print(dim(sce_SEG))
  # summary_SEG <- data.frame(mean_seg = colMeans(logcounts(sce_SEG)), row.names=colnames(sce_SEG))
  # # head(summary_SEG)
  # sce_SEG$mean_seg <- summary_SEG[colnames(sce_SEG),"mean_seg"]
  meta_data <- as.data.frame(colData(sce_SEG))
  meta_data$mean_seg <- colMeans(logcounts(sce_SEG))
  meta_data$batch_info <- as.factor(meta_data$batch_info)
  meta_data$library_label <- gsub("SCRNA10X_SA_CHIP","",meta_data$library_label)
  
  meta_data$group_info <- ifelse(meta_data$Grouping=="Metastasis","M",
                           ifelse(meta_data$Grouping=="Primary","P",""))
  meta_data$st_info <- ifelse(meta_data$Site_origin=="Axillary","Ax",
                              ifelse(meta_data$Site_origin=="Primary","Pr",
                                     ifelse(meta_data$Site_origin=="Supraspinal","Su",
                                            ifelse(meta_data$Site_origin=="Ventral_spinal","Ve",""))))
  meta_data$library_label <- paste0(meta_data$library_label,'_',meta_data$group_info,'_',meta_data$st_info)
  p <- plot_variation_function(meta_data, xstring="library_label", ystring="mean_seg", 
                               plottype="batch_info", plottitle,
                               xlabel='', ylabel="Mean exp of SE genes")
  if(!is.null(yl)){
    p <- p + ylim(yl[1],yl[2])
  }
  
  return(p)
  
}


plot_stable_exp_genes <- function(){
  input_dir = '/home/htran/storage/datasets/metastasis_results/rnaseq_SA919/'
  raw_sce <- readRDS(paste0(input_dir,'normalized/SA919_raw.rds'))
  normalized_sce <- readRDS(paste0(input_dir,'normalized/SA919_normalized_output.rds'))
  logcounts(raw_sce) <- log2(counts(raw_sce)+1)
  assayNames(normalized_sce)
  assayNames(raw_sce)
  dim(raw_sce)
  raw_sce <- raw_sce[rownames(normalized_sce),]
  
  assayNames(raw_sce)
  norm_df <- logcounts(sce2)
  dim(norm_df)
  norm_df <- logcounts(normalized_sce)
  dim(normalized_sce)
  
  var_genes <- setdiff(rownames(norm_df),obs_seg$gene_ens)
  length(t)
  classified_gene1 <- data.frame(gene_ens=c(obs_seg$gene_ens,var_genes),
                                 gene_type=c(rep('stable',length(obs_seg$gene_ens)),rep('variation',length(var_genes))),
                                 row.names = c(obs_seg$gene_ens,var_genes))
  
  ref_genes <- segList_ensemblGeneID$human$human_scSEG
  ref_genes <- intersect(rownames(normalized_sce), ref_genes)
  var_genes2 <- setdiff(rownames(norm_df),ref_genes)
  length(var_genes2)
  classified_gene2 <- data.frame(gene_ens=c(ref_genes, var_genes2),
                                 gene_type=c(rep('stable',length(ref_genes)),rep('variation',length(var_genes2))),
                                 row.names = c(ref_genes, var_genes2))
  
  dim(classified_gene2)
  
  p1 <- inspect_rawdata_std_mean(norm_df, classified_gene1)
  p1
  p2 <- inspect_rawdata_std_mean(norm_df, classified_gene2)
  p2
  png(paste0(input_dir,"batch_effect/scSEG_std_mean_plt_rawdata_33538g_504seg.png"), height = 2*500, width=2*700,res = 2*72)
  print(p1)
  dev.off()
  
  png(paste0(input_dir,"batch_effect/scSEG_std_mean_plt_rawdata_33538g_1024refgenes.png"), height = 2*500, width=2*700,res = 2*72)
  print(p2)
  dev.off()
  
  raw_sce$library_label[1:3]
  raw_sce$Grouping[1:2]
  summary(as.factor(raw_sce$Site_origin))
  pr <- get_plt_SEG(raw_sce, obs_seg$gene_ens,plottitle='Batch Effect Raw Data',yl=c(0,1))
  pr
  png(paste0(input_dir,"batch_effect/raw_data_scSEG.png"), height = 2*520, width=2*600,res = 2*72)
  print(pr)
  dev.off()
  pn <- get_plt_SEG(normalized_sce, obs_seg$gene_ens,plottitle='Batch Effect Normalized Data',yl=c(0,1)) 
  pn
  png(paste0(input_dir,"batch_effect/normalized_scSEG.png"), height = 2*520, width=2*600,res = 2*72)
  print(pn)
  dev.off()
  pn1 <- get_plt_SEG(normalized_sce, obs_seg$gene_ens,plottitle='Batch Effect Normalized Data',yl=NULL) 
  pn1
  png(paste0(input_dir,"batch_effect/normalized_scSEG2.png"), height = 2*520, width=2*600,res = 2*72)
  print(pn1)
  dev.off()
  ptotal <- cowplot::plot_grid(pr, pn, align = "h", ncol = 2)
  png(paste0(input_dir,"batch_effect/summary_scSEG.png"), height = 2*520, width=2*1200,res = 2*72)
  print(ptotal)
  dev.off()
  
  
  t = logcounts(raw_sce)
  t[1:3,1:3]
  plot_tsne_sce(raw_sce, normalized_sce, input_dir)
}

plot_tsne_sce <- function(raw_sce, normalized_sce, input_dir){
  raw_sce <- runPCA(raw_sce, ncomponents = 20, exprs_values = "logcounts")
  raw_sce <- runTSNE(raw_sce, dimred = "PCA", n_dimred = 20, ncomponents = 2)
  raw_sce <- runUMAP(raw_sce, dimred = "PCA", n_dimred = 20, ncomponents = 2, n_neighbors = 30)
  
  normalized_sce <- runPCA(normalized_sce, ncomponents = 20, exprs_values = "logcounts")
  normalized_sce <- runTSNE(normalized_sce, dimred = "PCA", n_dimred = 20, ncomponents = 2)
  normalized_sce <- runUMAP(normalized_sce, dimred = "PCA", n_dimred = 20, ncomponents = 2, n_neighbors = 30)
  
  pr_pca <- plotReducedDim(raw_sce, dimred="PCA", colour_by="batch_info")
  pn_pca <- plotReducedDim(normalized_sce, dimred="PCA", colour_by="batch_info")
  
  pr_umap <- plotReducedDim(raw_sce, dimred="UMAP", colour_by="batch_info")
  pn_umap <- plotReducedDim(normalized_sce, dimred="UMAP", colour_by="batch_info")
  pr_tsne <- plotReducedDim(raw_sce, dimred="TSNE", colour_by="batch_info")
  pn_tsne <- plotReducedDim(normalized_sce, dimred="TSNE", colour_by="batch_info")
  ppca <- cowplot::plot_grid(pr_pca, pn_pca, align = "h", ncol = 2)
  ptsne <- cowplot::plot_grid(pr_tsne, pn_tsne, align = "h", ncol = 2)
  pt <- cowplot::plot_grid(pr_pca, pr_tsne, pr_umap, pn_pca, pn_tsne, pn_umap, align = "hv", ncol = 3)
  png(paste0(input_dir,"batch_effect/reduction_plots.png"), height = 2*500, width=2*900,res = 2*72)
  print(pt)
  dev.off()
  # png(paste0(input_dir,"batch_effect/PCA.png"), height = 2*250, width=2*620,res = 2*72)
  # print(ppca)
  # dev.off()
  # png(paste0(input_dir,"batch_effect/TNSE.png"), height = 2*250, width=2*620,res = 2*72)
  # print(ptsne)
  # dev.off()
}


