# https://github.com/kstreet13/bioc2020trajectories/blob/d0a386f8c1411ebc40e6edc505531926c46f354e/vignettes/workshopTrajectories.Rmd#L230
# https://statomics.github.io/tradeSeq/articles/tradeSeq.html
# https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html#trajectory-inference
# https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html#trajectory-inference

# install.packages("BiocManager")
# BiocManager::install('SingleCellExperiment')
# BiocManager::install('dplyr')
# BiocManager::install('mclust')
# BiocManager::install('igraph')
# BiocManager::install('matrixStats')
# BiocManager::install('slingshot')
# BiocManager::install('stringr')
# BiocManager::install('tidyverse')
# BiocManager::install('RColorBrewer')
# BiocManager::install('scater')
# BiocManager::install('Seurat')

# BiocManager::install('rgl')
# BiocManager::install('tradeSeq')

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot, quietly = TRUE)
  library(mclust, quietly = TRUE)
  library(grDevices)
  library(RColorBrewer)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# Get sce file, and prepare data, and save into folder. 



 
head(cellData(crv1))
rownames(cellData(crv1))[1:3]
dim(pca_df)
length()
cells_use <- intersect(pca_df$cell_id, rownames(cellData(crv1)))
dim(crv1)
length(cells_use)
rownames(crv1)[1]
crv1 <- crv1[cells_use,]
pca_df <- data.table::fread('/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/hvg_umaps/SA609_norm_pca_10000.csv') %>% as.data.frame()
pca_df <- data.table::fread('/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/hvg_umaps/SA535_norm_pca_10000.csv') %>% as.data.frame()
rownames(pca_df) <- pca_df$cell_id
pca_df$cell_id[1:3]
pcas <- as.matrix(pca_df[rownames(cellData(crv1)),c('PC_1','PC_2')])
dim(pcas)
rownames(cellData(crv1))[1]
ls_lineages <- get_lineages(crv1, crv_umap_embed=NULL, pcas)
sce <- readRDS()
lg_df <- ls_lineages$lineages
umaps <- reducedDims(sce)[['UMAP']] %>% as.data.frame()
umaps <- as.data.frame(pcas)
rownames(umaps)[1]
umaps$cluster_label <- sce[,rownames(umaps)]$cluster_label
umaps$cluster_label <- as.factor(umaps$cluster_label)
class(umaps)
umaps$clone <- sce$clone
lg_df1 <- ls_lineages$unique_curves
length(unique(umaps$cluster_label))

umaps$treatmentSt <- colData(sce)[rownames(umaps),'treat']

# Draw clusters
res <- set_color_clusters(umaps$cluster_label)
lg <- plot_colors(res$clone_palette, output_dir, legend_label='clusters')
cols_use <- res$clone_palette
pc <- plot_curves(umaps, lg_df1, ls_lineages$curves, cols_use, output_dir, col_plt='cluster_label', 
                  x_plt='UMAP_1', y_plt='UMAP_2')
df_annot <- ls_lineages$curves
curves <- lg_df1
pc <- plot_curves(umaps, lg_df1, ls_lineages$curves, cols_use, output_dir, col_plt='cluster_label', 
                  x_plt='PC_1', y_plt='PC_2')

res2 <- set_color_treatmentSt(colData(sce)[,'treat'], datatag, cols=NULL)
lg2 <- plot_colors(res2$clone_palette, save_dir, legend_label='treatmentSt')
cols_use <- res2$clone_palette
pt <- plot_curves(umaps, lg_df1, ls_lineages$curves, cols_use, output_dir, col_plt='treatmentSt', 
                  x_plt='PC_1', y_plt='PC_2')

ptotal <- cowplot::plot_grid(pc, pt, ncol=2)
png(paste0(output_dir,"slingshot_output_",datatag,".png"), height = 2*500, width=2*900,res = 2*72)
print(ptotal)
dev.off()


# rd: data frame reduction dimensions
# curves: data frame
# rd <- umaps
# curves <- lg_df
# sds <- crv_umap_embed
# lid <- 'Lineage7'
# pt <- slingPseudotime(sds) %>% as.data.frame()













# datatag <- 'SA535'
# input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
# save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')
# sce_fn <- paste0(save_dir,datatag,'_raw.rds')
# sce <- prepare_data(sce_fn, datatag, save_dir, save_data=T, umaps=NULL, pcs=NULL, compute_PCA=T,compute_UMAP=T)
# pseudos <- get_slingshot_pseudotime(sce, save_dir, datatag, start_cls=NULL, cl=NULL, rd_use='PCA')
# "/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/"
# norm_sce <- readRDS(paste0(save_dir,'SA535_sctransform_normalized.rds'))
# dim(norm_sce)
# assayNames(norm_sce)
# sce <- norm_sce
# 
# ssd <- readRDS("/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/slingshot_pseudotime_SA535_start_8_UMAP_crv.rds")
# class(ssd)
# 
# pseudotime1 <- slingPseudotime(ssd, na = FALSE)
# dim(pseudotime1)
# View(head(pseudotime1))
# cellWeights <- slingCurveWeights(ssd)
# dim(cellWeights)
# View(head(cellWeights))
# ssd
# t <- SlingshotDataSet(ssd)
# class(ssd)
# plot(ssd, type = 'lineages', linI)
# 
# 
# slingLineages(ssd)
# getCurves(ssd)

# datatag <- 'SA609'
# save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')
# sce_fn <- paste0(save_dir,datatag,'_raw.rds')
# prepare_data(sce_fn, datatag, save_dir, save_data=T, umaps=NULL, pcs=NULL, compute_PCA=T,compute_UMAP=T)
# pseudos <- get_slingshot_pseudotime(sce, save_dir, datatag, start_cls=NULL, cl=NULL, rd_use='PCA')

# norm_mtx <- data.table::fread(paste0(save_dir,'sctransform_SA535.csv.gz'))
# dim(norm_mtx)
# genes_id <- data.table::fread(paste0(save_dir,'genes_ids_3000hvg.csv.gz'))
# dim(genes_id)
# rownames(norm_mtx) <- genes_id$gene_id
# colnames(norm_mtx)[1:4]
# norm_mtx$V1 <- NULL
# metacells <- data.table::fread(paste0(save_dir,'grouping_SA535_v2.csv.gz'))
# dim(metacells)
# colnames(metacells)
# metacells <- metacells %>%
#   as.data.frame()%>%
#   dplyr::select(-UMAP_1, -UMAP_2)
# 
# rownames(metacells) <- metacells$cell_id
# # View(head(metacells))
# sce <- SingleCellExperiment::SingleCellExperiment(assays=list(logcounts=as.matrix(norm_mtx)),
#                                                   colData = metacells[colnames(norm_mtx),])
# saveRDS(sce, sce_fn)
# dim(logcounts(sce))


# prepare_data <- function(sce_fn, datatag, save_dir, save_data=T, umaps=NULL, pcs=NULL, compute_PCA=F,compute_UMAP=F){
#   
#   if(!dir.exists(save_dir)){
#     dir.create(save_dir)
#   }
#   sce <- readRDS(sce_fn)
#   print(dim(sce))
#   if(compute_PCA & is.null(pcs)){
#     # pca_out <- prcomp(t(log1p(logcounts(sce))), scale. = FALSE)  
#     pca_out <- prcomp(t(logcounts(sce)), scale. = FALSE)
#     # rd <- pca$x[,1:2]
#     # dim(pca_out$x)
#     pcs <- pca_out$x[,1:15]
#   }
#   if(compute_UMAP & is.null(umaps)){
#     library(uwot)
#     rd <- uwot::umap(t(logcounts(sce)))
#     colnames(rd) <- c('UMAP_1', 'UMAP_2')
#   }  
#   # pcs <- pca_df[,paste0('PC_',rep(1:20,1))]
#   # umaps <- umap_df[,c('UMAP_1','UMAP_2')]
#   
#   ex_genes <- var_genes[!var_genes %in% rownames(sce)] # problem with _ and - 
#   # [1] "SOD2-ENSG00000112096"   
#   # [2] "POLR2J3-ENSG00000285437"
#   # [3] "TBCE-ENSG00000285053"   
#   # [4] "MATR3-ENSG00000015479"
#   genes_use <- c(intersect(var_genes, rownames(sce)),gsub('-','_',ex_genes))
#   sce <- sce[genes_use,]  
#   dim(sce)
#   reducedDims(sce) <- SimpleList(PCA = as.matrix(pcs), UMAP = as.matrix(umaps))
#   saveRDS(sce, paste0(save_dir, 'sce_raw_3000hvg.rds'))
#   cl <- mclust::Mclust(pcs)$classification
#   # umap_df <- data.table::fread(paste0(save_dir, 'SA535_norm_umap.csv')) %>% as.data.frame()
#   # dim(umap_df)
#   # head(umap_df)
#   # print(unique(cl))
#   # colnames(colData(sce))
#   colData(sce)$cluster_label <- cl
#   if(save_data){
#     saveRDS(sce, paste0(save_dir, 'slingshot_raw_',datatag,'.rds'))
#   }
#   sce <- readRDS(paste0(save_dir, 'sce_raw_3000hvg.rds'))
#   dim(sce)
#   return(sce)
# }


# get_slingshot_pseudotime <- function(sce, save_dir, datatag, 
#                                      start_cls=NULL, cl=NULL, rd_use='PCA'){
#   rd_use='UMAP'
#   start_cls <- '8'
#   rd <- reducedDims(sce)[[rd_use]]
#   dim(rd)
#   class(rd)
#   if(is.null(cl)){
#     cl <- colData(sce)[,'cluster_label'] # Gaussian mixture model clustering - optimize nb clusters  
#   }
#   
#   # head(rd)
#   print(unique(cl))
#   lin1 <- getLineages(rd, cl, start.clus = start_cls)
#   crv1 <- getCurves(lin1, approx_points = 150)
#   class(crv1)
#   cls <- unique(cl)
#   col <- brewer.pal(9,"Set1")
#   col <- col[1:length(cls)]
#   names(col) <- paste0('Cls_',cls)
#   
#   
#   
#   if(!is.null(start_cls)){
#     start_cls <- paste0('_start_',start_cls)
#     print(start_cls)
#     save_dir
#     
#   }
#   
#   saveRDS(crv1, paste0(save_dir, "slingshot_pseudotime_",datatag,start_cls,'_',rd_use,"_crv.rds"))
#   
#   # pt <- slingshot(rd, cl)
#   # pt
#   # sce <- slingshot(sce, clusterLabels = 'cluster_label',
#   #                  reducedDim = 'UMAP', approx_points = 100)
#   plot(rd, pch=16, asp = 1,
#        col = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[cl])
#   png(paste0(save_dir,"slingshot_",datatag,start_cls,'_',rd_use,"_clusters.png"), height = 2*600, width=2*500,res = 2*72)
#   plot(rd, pch=16, asp = 1, col = c(brewer.pal(9,"Set1")[cl]))
#   dev.off()     
#   
#   lines(SlingshotDataSet(pt), type = 'l', lwd=2, col='black')
#   length(unique(cl))
#   png(paste0(save_dir,"slingshot_",datatag,start_cls,'_',rd_use,"_lineage.png"), height = 2*400, width=2*500,res = 2*72)
#   plot(rd[,1:2], col = cols[sce$treat], asp = 1, pch = 16)
#   # lines(SlingshotDataSet(pt), type = 'lineages', lwd = 3, col = 'black')
#   lines(SlingshotDataSet(crv1), type = 'lineages', lwd = 3, col = 'black')
#   dev.off()
#   
#   png(paste0(save_dir,"slingshot_",datatag,start_cls,'_',rd_use,"_lineage_treatment.png"), height = 2*500, width=2*600,res = 2*72)
#   plot(rd[,1:2], col = cols[sce$treat], asp = 1, pch = 16)
#   # lines(SlingshotDataSet(pt), type = 'lineages', lwd = 3, col = 'black')
#   lines(SlingshotDataSet(crv1), type = 'lineages', lwd = 3, col = 'black')
#   dev.off()
#   
#   
#   meta_cells <- as.data.frame(colData(sce))
#   t <- table(meta_cells$treat,meta_cells$cluster_label)
#   t <- as.data.frame(t)
#   t <- t[t>=100]
#   col
#   t <- t[t$Freq>=100,]
#   for(f in unique(t$Var2)){
#     tmp <- t[t$Var2==f,]
#     print(paste0('Cls ',f,': ',paste(tmp$Var1, collapse = ',')))
#   }
#   cls1 <- paste0('Cls_',cl)
#   png(paste0(save_dir,"slingshot_",datatag,start_cls,'_',rd_use,"_ligneage.png"), height = 2*500, width=2*600,res = 2*72)
#   plot(rd[,1:2], col = col[cls1], asp = 1, pch = 16)
#   lines(SlingshotDataSet(crv1), type = 'lineages', lwd = 3, col = 'black')
#   # lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')
#   dev.off()
#   t@lineages
#   head(crv1)
#   colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
#   plotcol <- colors[cut(crv1, breaks=100)]
#   png(paste0(save_dir,"slingshot_pseudotime_",datatag,start_cls,".png"), height = 2*400, width=2*500,res = 2*72)
#   plot(rd[,1:2], col = plotcol, pch=16, asp = 1)
#   lines(SlingshotDataSet(crv1), lwd=2, col='black')  # smooth curve
#   dev.off()
#   
#   colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
#   plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
#   # SlingshotDataSet(sce)
#   sum(is.na(sce$slingPseudotime_1))
#   plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
#   dev.off()
#   plot(rd, pch=16, asp = 1)
#   lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')  # smooth curve
#   dim(sce)
#   
#   pseudos <- as.data.frame(colData(sce))  
#   pseudos$cell_id <- rownames(pseudos)
#   print(head(pseudos))  
#   data.table::fwrite(paste0(save_dir, "slingshot_pseudotime_",datatag,start_cls,".csv"))
#   return(pseudos)
# }  
# 
