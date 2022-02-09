# 
# # BiocManager::install('MOCCA')
# library(MOCCA)
# 
# set.seed(12345)
# data(toy5)
# res <- mocca(toy5, R=10, K=2:5)
# print(res$objectiveVals)
# # plot kmeans result for MCA index against neuralgas result for MCA index
# plot(res$objectiveVals[1,], res$objectiveVals[5,], pch=NA,
#      xlab=rownames(res$objectiveVals)[1], ylab=rownames(res$objectiveVals)[5])
# text(res$objectiveVals[1,], res$objectiveVals[5,], labels=colnames(res$objectiveVals))
# 
# pareto<-analyzePareto(res$objectiveVals)
# print(pareto)
# moo.clusterid<-pareto$rank[1]
# print(res$objectiveVals[,moo.clusterid])##find objective values for the combination of different clustering techniques and clustering indices measures
# res$cluster$kmeans[[moo.clusterid]][[10]]

# TO DO: DE between treated vs untreated cells --> DE genes 
# Genes clustering 

base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
gene_type <- 'increase'
save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/gene_regression_',gene_type,'_v2/')
save_dir
datatag <- 'SA609'
dir.create(save_dir, recursive = T, showWarning=F)
obs_clones <- c('R')
# clone_aware <- FALSE
clone_aware <- TRUE
obs_treatment_st <- c('UT','UTT','UTTT','UTTTT','UTTTTT')
obs_clones_untreated <- c('H')
obs_untreated_st <- c('UU','UUU','UUUU','UUUUU')
# sce <- readRDS(paste0(base_dir,'rnaseq_v6/SA609-v6/SA609_sctransform_normalized.rds'))

corrected_fn <- paste0(base_dir, 'rnaseq_v6/normalization_evaluation/',datatag,'/', datatag, "_norm_batch_corrected_sce.rds")
sce <- readRDS(corrected_fn)
sce_backup <- sce
dim(sce)
meta_data <- colData(sce) %>% as.data.frame()
class(meta_data)
summary(as.factor(meta_data$clone))
cells_use_g1 <- meta_data %>% 
  dplyr::filter(treatmentSt %in% obs_treatment_st &
                  clone %in% obs_clones) %>%
  dplyr::pull(cell_id)

cells_use_g2 <- meta_data %>% 
  dplyr::filter(treatmentSt %in% obs_untreated_st &
                  clone %in% obs_clones_untreated) %>%
  dplyr::pull(cell_id)

print(length(cells_use_g1))
cells_use_g1 <- cells_use_g1[sample(1:length(cells_use_g1), length(cells_use_g2), replace=FALSE)]
print(length(cells_use_g2))
sce <- sce[,c(cells_use_g1,cells_use_g2)]
srt <- Seurat::as.Seurat(sce, counts = "counts", data = "logcounts") 
meta_data <- srt@meta.data

input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
genes_map <- read.csv(paste0(input_dir, "biodatabase/meta_genes.csv"), header=T, stringsAsFactors=F)
rownames(genes_map) <- genes_map$gene_ens
dim(genes_map)


if(length(cells_use_g1) < 100 & length(cells_use_g2) < 100){
  print("There are no cells which satisfy the input condition or small nb cells only")
} else{
  groups_use <- c(obs_clones, obs_clones_untreated)
  feature_use="clone"
 
  calculate_DE_analysis_v2(paste0(datatag,'_cls',obs_clones,'_cls',obs_clones_untreated),
                           meta_data, 
                           cells_use_g1, cells_use_g2, feature_use,
                           srt, genes_map, groups_use, 'wilcox', 
                           save_dir, input_dir,
                           pAdjustThrs=0.05, minLogFC=0.25,
                           nbtopup = 30, nbtopdown = 30, save_data=T, viz=F)
  
  
  
}

data_dir <- paste0(save_dir,'SA609_R_H/')
deg <- read.table(paste0(data_dir,'de_significant_genes.txt'), header=T, sep='\t', check.names=F, stringsAsFactors=F, row.names=1)
deg_top <- read.table(paste0(data_dir,'de_topup_30_topdown_30.txt'), header=T, sep='\t', check.names=F, stringsAsFactors=F, row.names=1)
dim(deg)
View(head(deg_top))
summary(deg$avg_log2FC)
View(head(deg))
deg <- deg %>%
  dplyr::filter(avg_log2FC>=1)

observed_genes <- deg %>%
  dplyr::filter(avg_log2FC>=0.5)%>%
  dplyr::pull(GENEID)

length(observed_genes)

observed_genes <- deg %>%
  dplyr::filter(avg_log2FC<=-0.5)%>%
  dplyr::pull(GENEID)

# fsce <- logcounts(sce[observed_genes,])
# dim(fsce)
# set.seed(285);##set 285 as seed value ##all are good
# res <- mocca(as.matrix(fsce), R = 10, K = 2:10, iter.max = 10, nstart = 10)
# print(res$objectiveVals)
# # plot kmeans result for MCA index against neuralgas result for MCA index
# plot(res$objectiveVals[1,], res$objectiveVals[5,], pch=NA,
#      xlab=rownames(res$objectiveVals)[1], ylab=rownames(res$objectiveVals)[5])
# text(res$objectiveVals[1,], res$objectiveVals[5,], labels=colnames(res$objectiveVals))
# 
# pareto<-analyzePareto(res$objectiveVals)
# print(pareto)
# moo.clusterid<-pareto$rank[1]
# print(res$objectiveVals[,moo.clusterid])##find objective values for the combination of different clustering techniques and clustering indices measures
# res$cluster$kmeans[[moo.clusterid]][[10]]
# 
# sce <- sce_backup
mtx <- logcounts(sce[observed_genes,]) %>% as.data.frame()
dim(mtx)
# mtx$ID <- rownames(mtx)
# mtx <- mtx %>%
#   select(ID, everything())

# View(mtx[1:3,1:3])
data.table::fwrite(mtx, paste0(save_dir,'upRx.csv'), quote=F, row.names = T)

?data.table::fwrite
t <- as.matrix(t(logcounts(sce[observed_genes,])))
pcs <-  scater::calculatePCA(t, ncomponents=20)
dim(pcs)
geneTree = flashClust::flashClust(dist(pcs), method="average")

library(factoextra)

library("dbscan")
data("moons")
plot(moons, pch=20)
dim(moons)
View(head(moons))
cl <- hdbscan(moons, minPts = 5)
cl
plot(moons, col=cl$cluster+1, pch=20)
cl$hc
plot(cl$hc, main="HDBSCAN* Hierarchy")
