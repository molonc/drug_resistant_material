
suppressPackageStartupMessages({
  library(monocle3)
  library(dplyr)
  library(ggplot2)
  # library(pheatmap)
  library(TSCAN)
})

base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA535-v6/'
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/monocle_v2/'

datatag <- 'SA535'
cds <- readRDS(paste0(output_dir,datatag,'_cds_Rx_processed_v2.rds'))
sce_raw <- readRDS(paste0(base_dir,"scTransform_SA535_CIS_Rx.rds"))
dim(sce_raw)
sum(colnames(cds) %in% colnames(sce_raw))
rownames(cds)[1]
rownames(sce_raw)[1]
# rownames(sce_raw) <- rowData(sce_raw)$ID
dim(cds)
cells_use <- colnames(cds)[colData(cds)$cluster_label %in% c("Cls_1","Cls_25","Cls_3")]
cells_use[1]
norm_mtx <- logcounts(sce_raw[rownames(cds),cells_use])
dim(norm_mtx)
norm_mtx <- as.matrix(norm_mtx)
lpsmclust <- exprmclust(norm_mtx, clusternum = 3:5)

saveRDS(lpsmclust,paste0(output_dir,'tscan_clusters.rds'))
# lpsmclust <- readRDS(paste0(output_dir,'tscan_clusters.rds'))
p1 <- plotmclust(lpsmclust, show_cell_names = F, cell_name_size = 0.3)
png(paste0(output_dir,"tscan_clusters_",datatag,".png"), height = 2*450, width=2*500,res = 2*72)
print(p1)
dev.off()

lpsorder <- TSCANorder(lpsmclust)
saveRDS(lpsorder,paste0(output_dir,'tscan_lpsorder.rds'))
# lpsorder
diffval <- difftest(norm_mtx,lpsorder)
sum(diffval$qval < 0.05)
#Selected differentially expressed genes under qvlue cutoff of 0.05
head(row.names(diffval)[diffval$qval < 0.05])
# STAT2expr <- log2(lpsdata["STAT2",]+1)
p2 <- singlegeneplot(norm_mtx, TSCANorder(lpsmclust,flip=TRUE,orderonly=FALSE))

png(paste0(output_dir,"tscan_pseudotime_",datatag,".png"), height = 2*450, width=2*500,res = 2*72)
print(p2)
dev.off()

data.table::fwrite(diffval,paste0(output_dir,'tscan_diffval.csv'))
print('Completed')
