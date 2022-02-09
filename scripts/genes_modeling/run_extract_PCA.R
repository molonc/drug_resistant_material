library(factoextra)
library(data.table)
library(dplyr)
library(dbscan)

save_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/gene_regression_increase_v2/"
mtx <- data.table::fread(paste0(save_dir,'upRx.csv')) %>% as.data.frame()

rownames(mtx) <- mtx$V1
mtx$V1 <- NULL
rownames(mtx)[1]
colnames(mtx)[1]

dim(mtx)

res.pca <- prcomp(as.matrix(mtx), scale = F)
print(head(res.pca))
saveRDS(res.pca, paste0(save_dir,'PCAs.rds'))

cl <- hdbscan(res.pca$x, minPts = 5)
print(cl)
saveRDS(cl, paste0(save_dir,'cl_hdbscan.rds'))
print(summary(as.factor(cl$cluster)))
# plot(res.pca[1:2,], col=cl$cluster+1, pch=20)
# cl$hc
# plot(cl$hc, main="HDBSCAN* Hierarchy")

cl <- readRDS(paste0(save_dir,'cl_hdbscan.rds'))
cl$cluster
gene_cluster <- data.frame(ens_gene=rownames(mtx),
                           gene=rowData(observed_sce)$Symbol,
                           cluster=paste0('Cls_',cl$cluster))
head(gene_cluster)
