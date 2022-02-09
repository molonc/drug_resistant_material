library(dplyr)
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/'
datatag <- 'SA609'
sce <- readRDS(paste0(base_dir,'rnaseq_v6/SA609-v6/SA609_sctransform_normalized.rds'))
print(dim(sce))
if(is.null(rowData(sce)$Symbol)){
  genes_symb_df <- read.csv(paste0(base_dir,'biodatabase/meta_genes.csv'), check.names = F, stringsAsFactors = F)
  dim(genes_symb_df)
  rownames(genes_symb_df) <- genes_symb_df$gene_ens
  rowData(sce)$Symbol <- genes_symb_df[rownames(sce),'gene_symb']
}
rownames(sce)[1]
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/'
plot_umap <- function(sce, exprs, metadata, output_dir, datatag, dims = 1:15){
  # assay(sce, exprs)
  exprs <- 'logcounts'
  srt <- Seurat::as.Seurat(sce, counts = "counts", data = exprs)  
  # srt <- Seurat::FindVariableFeatures(object = srt, verbose = T)
  # srt <- Seurat::ScaleData(object = srt, verbose = FALSE) # problematics
  srt[["RNA"]]@scale.data <- as.matrix(logcounts(sce))
  # class(srt[["RNA"]]@data)
  nfeatures_use <- 5000
  srt <- Seurat::FindVariableFeatures(object = srt, verbose = FALSE, 
                                      selection.method = "vst", nfeatures = nfeatures_use)
  var_genes <- Seurat::VariableFeatures(object = srt)
  print(paste0('Nb var genes: ',length(var_genes)))
  var_genes_df <- data.frame(var_gene=var_genes)
  # var_genes_df$var_gene[1]
  write.csv(var_genes_df, file=paste0(output_dir,'var_genes_',nfeatures_use,'.csv'), quote=F, row.names = F)
  srt <- Seurat::RunPCA(object = srt, verbose = FALSE, npcs = 50) #, features = var_genes
  dims = 1:50
  srt <- Seurat::FindNeighbors(srt, dims = dims)
  print("Run UMAP")
  # # srt <- RunTSNE(object = srt, verbose = FALSE)
  srt <- Seurat::RunUMAP(object = srt, dims = dims, verbose = FALSE) #umap.method = "umap-learn", metric='correlation'
  p21 <- Seurat::DimPlot(srt, reduction = "umap", group.by = 'treatmentSt')
  p21 <- Seurat::DimPlot(srt, reduction = "pca", group.by = 'treatmentSt')
  
  pca_df <- as.data.frame(Seurat::Embeddings(object = srt, reduction = "pca"))
  umap_df <- as.data.frame(Seurat::Embeddings(object = srt, reduction = "umap"))
  meta_info <- srt@meta.data
  meta_info$cell_id <- rownames(meta_info)
  meta_info <- meta_info %>%
    dplyr::select(cell_id, Barcode, library_id, sample, treatmentSt, timepoint)
  
  umap_df$cell_id <- rownames(umap_df)
  umap_df <- umap_df %>% inner_join(meta_info, by=c("cell_id"))
  umap_df$cell_id <- paste0(umap_df$library_id,'_',umap_df$Barcode)
  
  pca_df$cell_id <- rownames(pca_df)
  pca_df <- pca_df %>% inner_join(meta_info, by=c("cell_id"))
  pca_df$cell_id <- paste0(pca_df$library_id,'_',pca_df$Barcode)
  dim(umap_df)
  dim(pca_df)
  
  data.table::fwrite(pca_df, paste0(output_dir, datatag,'_',nfeatures_use, "_norm_pca.csv"))
  data.table::fwrite(umap_df, paste0(output_dir, datatag,'_',nfeatures_use, "_norm_umap.csv"))
  
}  
colnames(sce)[1]
pca_df$cell_id[1]
pca_df$scell_id <- paste0(pca_df$sample,'_',pca_df$Barcode)
umap_df$scell_id <- paste0(umap_df$sample,'_',umap_df$Barcode)
sum(colnames(sce)==pca_df$scell_id)
sum(colnames(sce)==umap_df$scell_id)
colnames(pca_df)
reducedDims(sce) <- SimpleList(PCA = as.matrix(pca_df[,paste0('PC_',rep(1:50,1))]), UMAP = as.matrix(umap_df[,paste0('UMAP_',rep(1:2,1))]))
rd <- SimpleList(PCA = as.matrix(pca_df[,paste0('PC_',rep(1:50,1))]), UMAP = as.matrix(umap_df[,paste0('UMAP_',rep(1:2,1))]))
saveRDS(rd, paste0(output_dir, datatag,'_',nfeatures_use,'_rd.rds'))

dim(sce)
var_genes_df$var_gene[1]
sce_backup <- sce
nfeatures_use <- 3000
var_genes_df <- data.table::fread(paste0(output_dir,'var_genes_',nfeatures_use,'.csv')) %>% as.data.frame()
sce <- sce[var_genes_df$var_gene,]
dim(sce)
saveRDS(sce, paste0(output_dir, 'sce_raw_',nfeatures_use,'hvg.rds'))
pcs <- rd$PCA
dim(pcs)
library(mclust)
cl <- mclust::Mclust(pcs)$classification
saveRDS(cl,paste0(output_dir, 'optimal_clusters.rds'))
colData(sce)$cluster_label <- cl
metacells <- data.frame(colData(sce))
table(metacells$cluster_label, metacells$treatmentSt)
colData(sce)$cluster_label <- cl

rd_use='UMAP'
start_cls <- '8'
rd <- reducedDims(sce)[[rd_use]]
dim(rd)
class(rd)
if(is.null(cl)){
  cl <- colData(sce)[,'cluster_label'] # Gaussian mixture model clustering - optimize nb clusters  
}

# head(rd)
print(unique(cl))
lin1 <- getLineages(rd, cl, start.clus = start_cls)
crv1 <- getCurves(lin1, approx_points = 150)


clone_palette_20 <- c(
  "#be5f72", "#d74058", "#dc4229", "#a6552c", "#df956f", "#e47a33",
  "#d49f34", "#836e2c", "#b2ad5a", "#92b539", "#4c7d38", "#4dc041",
  "#5dba7f", "#47b8c3", "#6280ca", "#7b57db", "#ce8bd1", "#934f94",
  "#cb48cb", "#d74391"
)
cls_ls <- unique(sce$cluster_label)
cols <- clone_palette_20[1:length(cls_ls)]
names(cols) <- cls_ls
save_dir <- output_dir
png(paste0(save_dir,"slingshot_",datatag,start_cls,'_',rd_use,"_lineage.png"), height = 2*400, width=2*500,res = 2*72)
plot(rd[,1:2], col = cols[sce$cluster_label], asp = 1, pch = 16)
# lines(SlingshotDataSet(pt), type = 'lineages', lwd = 3, col = 'black')
lines(SlingshotDataSet(crv1), type = 'lineages', lwd = 3, col = 'black')
dev.off()
