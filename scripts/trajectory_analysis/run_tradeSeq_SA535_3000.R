# tradeSeq 

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot, quietly = TRUE)
  # library(mclust, quietly = TRUE)
  library(tradeSeq)
  library(grDevices)
  library(RColorBrewer)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})



datatag <- 'SA535'
input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')

output_dir <- paste0(save_dir,'tradeSeq_3000/')
dir.create(output_dir)
nfeatures_use <- 3000
sce <- readRDS(paste0(save_dir, datatag,'_',nfeatures_use,'_rd_sce.rds'))
norm_sce <- readRDS(paste0(save_dir, datatag,'_',nfeatures_use,'_rd_sce.rds'))
dim(norm_sce)
dim(sce)

# clone_df1 <- clone_df %>%
#   mutate(cell_id=paste0(library_id,'_',Barcode))%>%
#   select(cell_id, clone)
# sce$cell_id[1]
# sce$cell_id <- colnames(sce)
# sum(clone_df1$cell_id %in% sce$cell_id)
# clone_df1$cell_id[1]
# sce@colData <- as.data.frame(colData(sce)) %>%
#   left_join(clone_df1, by=c('cell_id'))%>%
#   DataFrame(check.names = FALSE)
# unique(sce$clone)
# sce$clone <- ifelse(is.na(sce$clone),'unassigned',sce$clone)
# saveRDS(sce, paste0(save_dir, datatag,'_',nfeatures_use,'_rd_sce_clones.rds'))
# unique(norm_sce$cluster_label)
# metacells <- colData(norm_sce)
# dim(metacells)
# metacells <- metacells %>%
#   as.data.frame() %>%
#   dplyr::filter(!cluster_label %in% c(7, 10, 8, 1, 6))
# rownames(metacells)[1]

var_genes_df <- data.table::fread(paste0(save_dir, "SA535_3000_hvg_genes.csv")) %>% as.data.frame()
# View(head(var_genes_df))
rownames(norm_sce)[1]
var_genes_df$bc[1]
norm_sce <- norm_sce[var_genes_df$bc[1:1000],] # get results first
counts <- as.matrix(counts(norm_sce))
dim(counts)
# meta_genes <- data.frame(ens_gene_id=rownames(norm_sce), gene_symbol=rowData(norm_sce)$Symbol, stringsAsFactors=F)
# rownames(meta_genes) <- meta_genes$ens_gene_id

# View(as.matrix(counts(sce)[1:5,1:5]))
# sce1 <- sce[, sce$cluster_label==6]
# unique(sce1$treat)
# table(sce$treat, sce$cluster_label)



# pseudotime <- data.table::fread(paste0(output_dir, "slingshot_SA535_7_PCA_pseudotime.csv")) %>% as.data.frame()
# # View(head(pseudo))
# pseudotime <- pseudotime %>% tibble::column_to_rownames(var='cell_id')
# 
# cellWeights <- data.table::fread(paste0(output_dir, "slingshot_SA535_7_PCA_cellWeights.csv")) %>% as.data.frame()
# cellWeights <- cellWeights %>% tibble::column_to_rownames(var='cell_id')
# View(head(cellWeights))
crv <- readRDS(paste0(save_dir, "slingshot_pseudotime_SA535_7_PCA_crv.rds"))
set.seed(7)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)

ts_sce <- tradeSeq::fitGAM(counts = counts,
                           pseudotime = pseudotime, cellWeights = cellWeights,
                           nknots = 7, verbose = FALSE,
                           parallel=T,
                           BPPARAM = BiocParallel::MulticoreParam(workers = 5))
saveRDS(ts_sce, paste0(output_dir, "fitGAM_out.rds"))
ts_sce <- readRDS(paste0(output_dir, "fitGAM_out.rds"))


patternRes <- tradeSeq::patternTest(ts_sce, l2fc = 0.5)
# oPat <- order(patternRes$waldStat, decreasing = TRUE)
# print(head(rownames(patternRes)[oPat]))
print(dim(patternRes))
saveRDS(patternRes, paste0(output_dir, "patternRes_out.rds"))
patternRes <- readRDS(paste0(output_dir, "patternRes_out.rds"))

endRes <- tradeSeq::diffEndTest(ts_sce, pairwise=TRUE, l2fc = 0.5)
saveRDS(endRes, paste0(output_dir, "endRes_out.rds"))
endRes <- readRDS(paste0(output_dir, "endRes_out.rds"))

earlyDERes <- tradeSeq::earlyDETest(ts_sce, l2fc = 0.5, pairwise=TRUE)
saveRDS(earlyDERes, paste0(output_dir, "earlyDERes_out.rds"))
earlyDERes <- readRDS(paste0(output_dir, "earlyDERes_out.rds"))




# Load data
output_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/tradeSeq_3000/tradeSeq_SA535/"
rev_genes <- data.table::fread(paste0(output_dir,'reversibility_genes/patternTest_diffEndTest_SA535_l2_Rx_lg1_RxH.csv'))
end_drug_genes <- data.table::fread(paste0(output_dir,'Rx_vs_UnRx/SA535_Rx_vs_UnRx_endTest.csv'))
early_drug_genes <- data.table::fread(paste0(output_dir,'Rx_vs_UnRx/SA535_Rx_vs_UnRx_earlyTest.csv'))
rev_genes_end <- data.table::fread(paste0(output_dir,'reversibility_genes/patternTest_diffEndTest_SA535_l2_Rx_lg1_RxH.csv'))
rev_genes_early <- data.table::fread(paste0(output_dir,'reversibility_genes/patternTest_earlyDERes_SA535_l2_Rx_lg1_RxH.csv'))

thrs_stat <- 200
statSecondTest_thrs <- 40
dim(rev_genes)
dim(rev_genes_end)
dim(rev_genes_early)
summary(rev_genes$pattern)
rev_genes_end <- rev_genes_end %>%
  filter(pattern>thrs_stat & statSecondTest>statSecondTest_thrs)
rev_genes_early <- rev_genes_early %>%
  filter(pattern>thrs_stat & statSecondTest>statSecondTest_thrs)
dim(rev_genes_end)
dim(rev_genes_early)
rev_genes1 <- dplyr::bind_rows(rev_genes_end, rev_genes_early)%>%
              dplyr::arrange(desc(pattern))
rev_genes1 <- rev_genes %>%
  filter(pattern>thrs_stat & statSecondTest>statSecondTest_thrs)%>%
  # filter(statSecondTest>mean(statSecondTest))%>%
  dplyr::arrange(desc(pattern))
# rev_genes <- rev_genes[1:50,]
dim(rev_genes1)
View(rev_genes1)
data.table::fwrite(rev_genes1, paste0(output_dir,'reversibility_genes/patternTest_diffEndTest_SA535_l2_Rx_lg1_RxH_topgenes_early_end.csv'))

# sa535_genes <- rev_genes1$gene_symbol
# intersect(sa535_genes, sa609_genes)
# summary(end_drug_genes$pattern)
# summary(early_drug_genes$pattern)

dim(end_drug_genes)
dim(early_drug_genes)
end_drug_genes <- end_drug_genes %>%
  filter(pattern>thrs_stat)%>%
  dplyr::arrange(desc(pattern))
early_drug_genes <- early_drug_genes %>%
  filter(pattern>thrs_stat)%>%
  dplyr::arrange(desc(pattern))
dim(end_drug_genes)
dim(early_drug_genes)

dim(rev_genes)
genes_use <- union(end_drug_genes$ens_gene_id, early_drug_genes$ens_gene_id)

genes_use <- unique(c(end_drug_genes$ens_gene_id, early_drug_genes$ens_gene_id, rev_genes_end$ens_gene_id))

end_drug_genes$ens_gene_id[1]
early_drug_genes$ens_gene_id[1]
rev_genes_end$ens_gene_id[1]
length(genes_use)

# genes_use <- genes_use[!genes_use %in% rev_genes1$ens_gene_id]
length(genes_use)
genes_use[1]
colnames(colData(ts_sce))
colData(ts_sce)$slingshot[1:5]

genes_use <- obs_genes_df$genes_use
preprocess_mat <- as.matrix(counts(ts_sce[genes_use,]))

dim(rev_genes)
# preprocess_mat <- as.matrix(counts(ts_sce[rev_genes$ens_gene_id,]))
# preprocess_mat[1:3,1:3]

gene_module_df <- get_genes_modules(preprocess_mat=preprocess_mat, resolution=0.1)
summary(as.factor(gene_module_df$module))
data.table::fwrite(gene_module_df, paste0(output_dir,'reversibility_genes/gene_modules_',datatag,'_transient.csv'))
dim(gene_module_df)
data.table::fwrite(gene_module_df, paste0(output_dir,'Rx_vs_UnRx/gene_modules_',datatag,'.csv'))
dim(gene_module_df)
gene_module_df <- data.table::fread(paste0(output_dir,'Rx_vs_UnRx/gene_modules_',datatag,'.csv'))
View(head(gene_module_df))
# obs_genes_df <- data.frame(genes_use=c(gene_module_df$id, rev_genes$ens_gene_id),
#                            gene_type=c(paste0('Cls',gene_module_df$module),
#                                        rep('transient_genes',length(rev_genes$ens_gene_id))))
# 

obs_genes_df <- data.frame(genes_use=gene_module_df$id,
                           gene_type=paste0('cls_',gene_module_df$module), 
                           row.names=gene_module_df$id)
# unique(gene_module_df$module)
# rev <- gene_module_df %>%
#     dplyr::filter(module %in% c(1, 4))
# sum(rev$id %in% rev_genes$ens_gene_id)
# sum(rev$id %in% compare_df$ens_gene_id)
# dim(rev)

dim(rev_genes)
test <- rev_genes %>%
  dplyr::filter(ens_gene_id %in% rev$id)
summary(test$statSecondTest)
obs_genes_df <- data.frame(genes_use=gene_module_df$id,
                           gene_type=paste0('cls_',gene_module_df$module), 
                           row.names=gene_module_df$id)
dim(obs_genes_df)
obs_genes_df <- obs_genes_df %>%
            dplyr::filter(gene_type %in% c('cls_2','cls_3'))
obs_genes_df <- obs_genes_df %>%
  dplyr::filter(!gene_type %in% c('cls_2','cls_3'))

# Repressed
obs_genes_df <- obs_genes_df %>%
  dplyr::filter(gene_type %in% c('cls_3'))
obs_genes_df <- obs_genes_df %>%
  dplyr::filter(gene_type %in% c('cls_2'))
dim(obs_genes_df)
rownames(patternRes)[1]
patternRes$ens_gene_id <- rownames(patternRes)
rp <- patternRes %>%
  filter(ens_gene_id %in% obs_genes_df$genes_use)%>%
  dplyr::arrange(desc(waldStat))
head(patternRes)
head(rp)
obs_genes_df <- rp
# head(obs_genes_df)

## SA609
# repressed_cls <- paste0('Cls',c(6,1))
# activated_cls <- paste0('Cls',c(2,3,5,4,7))
repressed_cls <- paste0('cls_',c(3))
activated_cls <- paste0('cls_',c(2))

obs_genes_df$gene_type <- ifelse(obs_genes_df$gene_type %in% repressed_cls,'Repressed',
                                 ifelse(obs_genes_df$gene_type %in% activated_cls,'Activated','Transient'))

summary(as.factor(obs_genes_df$gene_type))
obs_genes_df <- obs_genes_df %>%
  dplyr::filter(gene_type=='Repressed')

t <- patternRes%>%
  dplyr::filter()
dim(obs_genes_df)
tag <- 'clusters_activated_repressed'
p <- viz_heatmap(ts_sce, obs_genes_df, output_dir, datatag, tag)

tag <- 'clusters_transient'
head(obs_genes_df)
rownames(obs_genes_df) <- obs_genes_df$genes_use
p1 <- viz_heatmap(ts_sce, obs_genes_df, output_dir, datatag, tag)
dim(obs_genes_treatment)
obs_genes_df$gene_type <- ifelse(obs_genes_df$gene_type %in% repressed_cls,'Repressed',
                                 ifelse(obs_genes_df$gene_type %in% activated_cls,'Activated','Transient'))
summary(as.factor(obs_genes_df$gene_type))
s <- obs_genes_df %>%
  group_by(gene_type)%>%
  summarise(nb_gene=n())%>%
  mutate(gene_type_desc=paste0(gene_type,'(',nb_gene,')'))%>%
  as.data.frame()
obs_genes_df <- obs_genes_df %>% inner_join(s, by='gene_type')
dim(obs_genes_df)
unique(obs_genes_df$gene_type)
obs_genes_df$gene_type <- obs_genes_df$gene_type_desc
# %>%
#   select(-gene_type)
# summary(as.factor(obs_genes_df$gene_type_desc))
head(meta_genes$gene_symbol)
obs_genes <- c('MYC','HIST1H4C','MKI67','ODAM')
obs_genes <- c('S100A6','SAA1')
ext_genes <- meta_genes %>%
  dplyr::filter(gene_symbol %in% obs_genes)
obs_genes <- ext_genes$ens_gene_id
# grep('MYC',meta_genes$gene_symbol, value = T)

obs_genes <- c('CD99','ODAM')
head(meta_genes)


input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/tradeSeq_3000/tradeSeq_SA535/'
genes_df <- readRDS(paste0(input_dir, 'clusters_activated_repressed.rds'))
dim(genes_df)
meta_genes <- genes_df@matrix_param$row_split
class(meta_genes)

class(genes_df)
genes_df
meta_genes <- meta_genes %>%
  select(gene_type)
meta_genes$ens_gene_id <- rownames(meta_genes)
data.table::fwrite(meta_genes, paste0(input_dir, 'genes_activated_repressed_SA535.csv'))
transient_genes_df <- readRDS(paste0(input_dir, 'clusters_transient.rds'))
dim(transient_genes_df)
meta_genes <- transient_genes_df@matrix_param$row_split
data.table::fwrite(meta_genes, paste0(input_dir, 'transient_genes_SA535.csv'))
dim(meta_genes)




trans_df <- data.table::fread(paste0(input_dir, 'transient_genes_SA535.csv')) %>% as.data.frame()
dim(trans_df)
unique(trans_df$gene_type)
trans_df$gene_type <- gsub('cls_','Trans',trans_df$gene_type)
# trans_df$gene_type <- paste0('Trans(',dim(trans_df)[1],')')
colnames(trans_df)

activated_repressed_df <- data.table::fread(paste0(input_dir, 'genes_activated_repressed_SA535.csv')) %>% as.data.frame()
dim(activated_repressed_df)
unique(activated_repressed_df$gene_type)
colnames(activated_repressed_df)

total_genes <- dplyr::bind_rows(activated_repressed_df, trans_df)
total_genes <- total_genes %>% left_join(annots, by=c('ens_gene_id'))
data.table::fwrite(total_genes, paste0(output_dir, datatag,'_total_genes_modules_act_repr_trans.csv'))
total_genes <- total_genes %>%
  dplyr::rename(genes_use=ens_gene_id)
plttitle <- 'Activated Repressed Transient genes'
p <- viz_heatmap(ts_sce, total_genes, save_figs_dir, datatag, plttitle)
p


# Get random genes and plot the trajectories
patternRes$ens_gene_id <- rownames(patternRes)
activated_genes <- activated_repressed_df %>% 
  dplyr::filter(grepl('Activated',gene_type)) %>% 
  left_join(patternRes,by='ens_gene_id')%>%
  dplyr::arrange(desc(waldStat))%>%
  dplyr::pull(ens_gene_id)

gact_ls <- activated_genes[1:5]

repressed_genes <- activated_repressed_df %>% 
  dplyr::filter(grepl('Repressed',gene_type)) %>% 
  left_join(patternRes,by='ens_gene_id')%>%
  dplyr::arrange(desc(waldStat))%>%
  dplyr::pull(ens_gene_id)
gpr_ls <- repressed_genes[1:5]

trans_df <- trans_df %>% 
  left_join(patternRes,by='ens_gene_id') 

dim(trans_df)
unique(trans_df$gene_type)
trans_genes1 <- trans_df %>%
  dplyr::filter(grepl('Trans1',gene_type)) %>%
  dplyr::arrange(desc(waldStat))%>%
  dplyr::pull(ens_gene_id)

gm1_ls <- trans_genes1[1:3]

trans_genes2 <- trans_df %>%
  dplyr::filter(grepl('Trans2',gene_type)) %>%
  dplyr::arrange(desc(waldStat))%>%
  dplyr::pull(ens_gene_id)

gm2_ls <- trans_genes2[1:3]

trans_genes3 <- trans_df %>%
  dplyr::filter(grepl('Trans3',gene_type)) %>%
  dplyr::arrange(desc(waldStat))%>%
  dplyr::pull(ens_gene_id)

gm3_ls <- trans_genes3[1:3]

# plt_ls[[count]] <- p1
plt_ls <- list()
obs_genes_ls <- c(gact_ls, gpr_ls, gm1_ls, gm2_ls, gm3_ls)
colnames(total_genes)
names(obs_genes_ls) <- c(rep('activated', length(gact_ls)),
                         rep('repressed', length(gpr_ls)),
                         rep('transient1', length(gm1_ls)),
                         rep('transient2', length(gm2_ls)),
                         rep('transient3', length(gm3_ls)))
figs_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/figs_v3/"
total_genes <- as.data.frame(total_genes)
total_genes$genes_use


meta_genes <- total_genes %>%
  dplyr::filter(genes_use %in% obs_genes_ls)
library(tidyverse)
meta_genes <- meta_genes %>% 
  remove_rownames %>%
  dplyr::select(genes_use, gene_symbol) %>%
  tibble::column_to_rownames('gene_symbol')
meta_genes1 <- meta_genes
meta_genes1$gene_symbol <- rownames(meta_genes1)
rownames(meta_genes1) <- meta_genes1$genes_use
obs_genes_ls1 <- meta_genes1[obs_genes_ls,'gene_symbol']
names(obs_genes_ls1) <- names(obs_genes_ls)
viz_given_gene_exp_lineages(obs_genes_ls1[11:length(obs_genes_ls1)], meta_genes, ts_sce, save_figs_dir, datatag)




# Change the name of the modules
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/tradeSeq_3000/'
genes_df <- data.table::fread(paste0(save_dir, 'SA535_total_genes_modules_act_repr_trans.csv')) %>% as.data.frame()
dim(genes_df)
summary(as.factor(genes_df$gene_type))

meta_types <- data.frame(gene_type=unique(genes_df$gene_type),
                         gene_type_module=c('Module5','Module4','Module1',
                                            'Module3','Module2'))
genes_df <- genes_df %>% inner_join(meta_types, by='gene_type')
colnames(genes_df)
genes_df <- genes_df %>%
  dplyr::select(ens_gene_id, gene_symbol, gene_type_module, everything())
  

head(genes_df)

data.table::fwrite(genes_df, paste0(save_dir, 'SA535_total_genes_modules_act_repr_trans_08_Dec.csv.gz'))

genes_df <- data.table::fread(paste0(save_dir, 'SA535_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
dim(genes_df)
colnames(genes_df)
patternRes <- readRDS('/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/tradeSeq_3000/patternRes_out.rds')
dim(patternRes)
patternRes$ens_gene_id <- rownames(patternRes)
save_dir_pw <- paste0(save_dir,'pathway_results/')
dir.create(save_dir_pw)
datatag <- 'SA535'
patternRes1 <- annotables::grch38 %>%
  # dplyr::select(gene_symbol = symbol, chr, start, end) %>%
  dplyr::select(gene_symbol = symbol, ens_gene_id = ensgene) %>%
  inner_join(patternRes, by='ens_gene_id')

deg_df <- genes_df %>% inner_join(patternRes, by='ens_gene_id')
dim(deg_df)
gsea_out <- get_pathway_trajectory_genes(patternRes1, save_dir_pw, desc='total_genes_modules', base_name = datatag, pathway_name='hallmark')
gsea_out <- get_pathway_trajectory_genes(patternRes1, save_dir_pw, desc='total_genes_modules', base_name = datatag, pathway_name='kegg')
head(gsea_out)
gsea_out$pathway


## Visualize specific genes expr
datatag <- 'SA535'
figs_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/figs_v3/"
# 
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
total_genes <- data.table::fread(paste0(save_dir, datatag,'_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
dim(total_genes)
obs_genes_ls <- total_genes$ens_gene_id
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/'
ts_sce <- readRDS(paste0(input_dir,'tradeSeq_3000/', "fitGAM_out.rds"))
dim(ts_sce)
output_fn <- paste0(input_dir,'tradeSeq_3000/','sigf_gene_exp.csv.gz')
avg_gene_exp <- get_average_gene_exp_per_lineage(datatag, ts_sce, obs_genes_ls, output_fn, nEstimatedPoints=100, save_data=T)
# save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/')
obs_genes_ls <- c('SAA1','CENPW','SMC4','RRM2')
meta_genes <- total_genes %>%
  dplyr::filter(gene_symbol %in% obs_genes_ls)

library(tidyverse)
meta_genes <- meta_genes %>% 
  remove_rownames %>%
  dplyr::select(ens_gene_id, gene_symbol) %>%
  tibble::column_to_rownames('gene_symbol')

head(meta_genes)
viz_given_gene_exp_lineages(obs_genes_ls, meta_genes, avg_gene_exp, figs_dir, datatag)





## Heatmap average expression plot
## Extracting heatmap of average genes expression x lineages 
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
genes_df <- data.table::fread(paste0(save_dir, 'SA535_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
dim(genes_df)
colnames(genes_df)
# summary(as.factor(genes_df$gene_type_module))
genes_df <- genes_df %>%
  dplyr::select(-gene_type) %>%
  dplyr::rename(gene_type=gene_type_module)
summary(as.factor(genes_df$gene_type))
plttitle <- 'Activated Repressed Transient genes'
save_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/"
ts_sce <- readRDS(paste0(save_dir,'tradeSeq_3000/', "fitGAM_out.rds"))
dim(ts_sce)
phm <- viz_heatmap(ts_sce, genes_df, save_dir, datatag, plttitle)
output_dir <- save_dir
exp_mtx <- data.table::fread(paste0(output_dir,"mtx_hm.csv.gz")) %>% as.data.frame()
obs_genes_df <- data.table::fread(paste0(output_dir,"obs_genes_hm.csv.gz")) %>% as.data.frame()
obs_cells_df <- data.table::fread(paste0(output_dir,"obs_cells_hm.csv.gz")) %>% as.data.frame()
rownames(exp_mtx) <- exp_mtx$ens_gene_id
exp_mtx$ens_gene_id <- NULL
dim(exp_mtx)
dim(obs_genes_df)
dim(obs_cells_df)
input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/cis_trans_lineages/')
cistrans_anno <- data.table::fread(paste0(input_dir, datatag,'_genes_cis_trans_lineages.csv')) %>% as.data.frame()
dim(cistrans_anno)
# head(cistrans_anno)
rownames(cistrans_anno) <- NULL
cistrans_anno$gene_type <- paste0("in ",cistrans_anno$gene_type)
meta_clone_lg <- data.table::fread(paste0(input_dir,datatag,'_meta_clone_lineages.csv')) %>% as.data.frame()
meta_clone_lg <- meta_clone_lg[!duplicated(meta_clone_lg$lineage),]

# sum(obs_genes_df$ens_gene_id==rownames(exp_mtx))
unique(obs_genes_df$gene_type)
gap <- viz_genes_exp_lineages_cistrans_anno_hm(as.matrix(exp_mtx), obs_genes_df, obs_cells_df,
                                               cistrans_anno, meta_clone_lg,
                                               output_dir, plttitle)

