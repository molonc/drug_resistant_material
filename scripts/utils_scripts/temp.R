library(Seurat)

input_dir <- '/home/htran/storage/rnaseq_datasets/hakwoo_metastasis_RNAseq/SA535_human/'
lid <- 'SCRNA10X_SA_CHIP0250_001'
sce <- readRDS(paste0(input_dir, lid,'/',lid,'.rdata'))
input_fn <- "/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/SCRNA10X_SA_CHIP0250_001_filtered.rds"
sce <- readRDS(input_fn)
dim(sce)
srt <- Seurat::as.Seurat(sce, counts = "counts", data=NULL)
res = 0.3
dims = 1:15
nb_hvg = 3000
srt <- Seurat::FindVariableFeatures(srt, selection.method = "vst", nfeatures = nb_hvg)
srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 100000)
srt <- ScaleData(srt, features = rownames(srt))
srt <- RunPCA(object = srt, verbose = FALSE)
srt <- FindNeighbors(srt, dims = dims)
srt <- FindClusters(srt, resolution = res)
srt <- RunUMAP(object = srt, dims = dims, verbose = FALSE) #umap.method = "umap-learn", metric='correlation'
# 
p1 <- DimPlot(srt, reduction = "umap")
p1
cls <- srt@meta.data$seurat_clusters
t[1:40]
cls <- srt$seurat_clusters
cls[1:30]
summary(as.factor(cls))
ggsave(
  plot = p1,
  height = 6,
  width = 8,
  dpi = 250,
  # To Do: add height, width here
  filename = paste0("/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/filtered/",lid,"_umap.png"),
  bg = "transparent"
)  
