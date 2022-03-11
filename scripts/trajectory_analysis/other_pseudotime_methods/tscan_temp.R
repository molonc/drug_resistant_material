library(TSCAN)
?plotmclust
?difftest
lpsmclust <- readRDS(paste0(output_dir,'tscan_clusters.rds'))
p1 <- plotmclust(lpsmclust, show_cell_names = F, cell_name_size = 0.3)
png(paste0(output_dir,"tscan_clusters_",datatag,".png"), height = 2*450, width=2*500,res = 2*72)
print(p1)
dev.off()

lpsmclust2 <- lpsmclust
length(lpsmclust2$clusterid)
cells_rd <- lpsmclust2$pcareduceres
dim(cells_rd)
t <- reducedDims(cds[,rownames(cells_rd)])[['UMAP']]
rownames(t)[1]
class(t)
lpsmclust2$pcareduceres <- reducedDims(cds[,rownames(cells_rd)])[['UMAP']]
class(cells_rd)
colnames(cells_rd)
rownames(cells_rd)[1]
ts <- colData(cds[,rownames(cells_rd)])$treatmentSt
ts[1]
lpsmclust2$clusterid <- ts
length(lpsmclust2$clusterid)
lpsmclust$clusterid[1:3]
names(lpsmclust2$clusterid) <- rownames(cells_rd)
lpsmclust3 <- lpsmclust2
lpsmclust3$clusterid <- lpsmclust$clusterid
names(lpsmclust3$clusterid)[1]
p13 <- plotmclust(lpsmclust3, show_cell_names = F, cell_name_size = 0.3, show_tree=F)

p11 <- plotmclust(lpsmclust2, show_cell_names = F, cell_name_size = 0.3, show_tree=F)
p11 <- plotmclust(lpsmclust2, show_cell_names = F, cell_name_size = 0.3)
png(paste0(output_dir,"tscan_treatmentst_",datatag,".png"), height = 2*450, width=2*500,res = 2*72)
print(p11)
dev.off()

p_tscan <- cowplot::plot_grid(p1, p11, ncol = 2)
p_tscan <- cowplot::plot_grid(p13,p11, ncol = 2)
png(paste0(output_dir,"tscan_trajectory_",datatag,"_umap.png"), height = 2*450, width=2*1000,res = 2*72)
print(p_tscan)
dev.off()



lpsorder <- readRDS(paste0(output_dir,'tscan_lpsorder.rds'))
p2 <- singlegeneplot(norm_mtx[100,], TSCANorder(lpsmclust,flip=TRUE,orderonly=FALSE))
class(norm_mtx)
png(paste0(output_dir,"tscan_pseudotime_",datatag,".png"), height = 2*450, width=2*500,res = 2*72)
print(p2)
dev.off()