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
  library(BiocParallel)
})



datatag <- 'SA609'
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/')
nfeatures_use <- 3000
output_dir <- paste0(save_dir,'tradeSeq_3000/')
# dir.create(output_dir)
norm_sce <- readRDS(paste0(save_dir, datatag,'_',nfeatures_use,'_rd_sce.rds'))
dim(norm_sce)

var_genes_df <- data.table::fread(paste0(save_dir, datatag,'_',nfeatures_use,"_hvg_genes.csv")) %>% as.data.frame()
# View(head(var_genes_df))
rownames(norm_sce)[1]
var_genes_df$bc[1]
norm_sce <- norm_sce[var_genes_df$bc[1:1000],] # get results first
dim(norm_sce)
# View(as.matrix(counts(sce)[1:5,1:5]))
# sce1 <- sce[, sce$cluster_label==6]
# unique(sce1$treat)
# table(sce$treat, sce$cluster_label)

counts = as.matrix(counts(norm_sce))

# meta_genes <- data.frame(ens_gene_id=rownames(norm_sce), gene_symbol=rowData(norm_sce)$Symbol, stringsAsFactors=F)
# rownames(meta_genes) <- meta_genes$ens_gene_id

# pseudotime <- data.table::fread(paste0(output_dir, "slingshot_SA535_7_PCA_pseudotime.csv")) %>% as.data.frame()
# # View(head(pseudo))
# pseudotime <- pseudotime %>% tibble::column_to_rownames(var='cell_id')
# 
# cellWeights <- data.table::fread(paste0(output_dir, "slingshot_SA535_7_PCA_cellWeights.csv")) %>% as.data.frame()
# cellWeights <- cellWeights %>% tibble::column_to_rownames(var='cell_id')
# View(head(cellWeights))
crv <- readRDS(paste0(save_dir, "slingshot_pseudotime_SA609_10_PCA_crv.rds"))
set.seed(7)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
ts_sce <- tradeSeq::fitGAM(counts = as.matrix(counts(norm_sce)),
                           pseudotime = pseudotime, cellWeights = cellWeights,
                           nknots = 8, verbose = FALSE,
                           parallel=T,
                           BPPARAM = BiocParallel::MulticoreParam(workers = 5))
saveRDS(ts_sce, paste0(output_dir, "fitGAM_out.rds"))
ts_sce <- readRDS(paste0(save_dir,'tradeSeq_3000/', "fitGAM_out.rds"))

patternRes <- tradeSeq::patternTest(ts_sce, l2fc = 0.5)
# oPat <- order(patternRes$waldStat, decreasing = TRUE)
# print(head(rownames(patternRes)[oPat]))
print(dim(patternRes))
saveRDS(patternRes, paste0(output_dir, "patternRes_out.rds"))
patternRes <- readRDS(paste0(save_dir,'tradeSeq_3000/', "patternRes_out.rds"))
# dim(patternRes)
colnames(patternRes)

endRes <- tradeSeq::diffEndTest(ts_sce, pairwise=TRUE, l2fc = 0.5)
saveRDS(endRes, paste0(output_dir, "endRes_out.rds"))


meta_genes <- data.frame(genes_use=patternRes$ens_gene_id,
                         gene_symbol=patternRes$gene_symbol, row.names = patternRes$gene_symbol)
dim(meta_genes)
data.table::fwrite(meta_genes, paste0(output_dir,'meta_genes.csv'))

# assoRes <- associationTest(ts_sce, l2fc = 0.5)
# saveRDS(assoRes, paste0(output_dir, "assoRes_out.rds"))
# assoRes <- readRDS(paste0(output_dir, "assoRes_out.rds"))


# startRes <- tradeSeq::startVsEndTest(ts_sce, lineages=TRUE, l2fc = 0.5)
# saveRDS(startRes, paste0(output_dir, "startRes_out.rds"))
# print(head(startRes))

# startRes <- readRDS(paste0(output_dir, "startRes_out.rds"))


# endRes <- readRDS(paste0(output_dir, "endRes_out.rds"))

earlyDERes <- tradeSeq::earlyDETest(ts_sce, l2fc = 0.5, pairwise=TRUE)
# ?tradeSeq::earlyDETest
saveRDS(earlyDERes, paste0(output_dir, "earlyDERes_out.rds"))
class(earlyDERes)
head(earlyDERes)
dim(earlyDERes1)

rownames(ts_sce)[1]
### based on mean smoother






# Load data
datatag <- 'SA609'
output_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/"
rev_genes <- data.table::fread(paste0(output_dir,'reversibility_genes/patternTest_diffEndTest_SA609_l1_Rx_lg2_RxH.csv'))
end_drug_genes <- data.table::fread(paste0(output_dir,'Rx_vs_UnRx/SA609_Rx_vs_UnRx_endTest.csv'))
early_drug_genes <- data.table::fread(paste0(output_dir,'Rx_vs_UnRx/SA609_Rx_vs_UnRx_earlyTest.csv'))
thrs_stat <- 200
statSecondTest_thres <- 40

rev_genes <- rev_genes %>%
  filter(pattern>thrs_stat & statSecondTest>statSecondTest_thres)%>%
  # filter(statSecondTest>mean(statSecondTest))%>%
  dplyr::arrange(desc(pattern))
dim(rev_genes)
sum()
rev_genes <- rev_genes %>%
  filter(pattern>thrs_stat)%>%
  # filter(statSecondTest>mean(statSecondTest))%>%
  dplyr::arrange(desc(pattern))
# rev_genes <- rev_genes[1:30,]  # get top 30 genes

data.table::fwrite(rev_genes, paste0(output_dir,'reversibility_genes/patternTest_diffEndTest_SA609_l1_Rx_lg2_RxH_top50.csv'))
rev_genes <- data.table::fread(paste0(output_dir,'reversibility_genes/patternTest_diffEndTest_SA609_l1_Rx_lg2_RxH_top50.csv'))
dim(rev_genes)
View(head(rev_genes))
summary(end_drug_genes$waldStat_global_endTest)
summary(early_drug_genes$pattern)


end_drug_genes <- end_drug_genes %>%
  filter(pattern>thrs_stat & waldStat_global_endTest>statSecondTest_thres)%>%
  dplyr::arrange(desc(pattern))
early_drug_genes <- early_drug_genes %>%
  filter(pattern>thrs_stat & waldStat_global_earlyTest>statSecondTest_thres)%>%
  dplyr::arrange(desc(pattern))
dim(end_drug_genes)
dim(early_drug_genes)
summary(early_drug_genes$pvalue_global_earlyTest)

genes_use <- union(end_drug_genes$ens_gene_id, early_drug_genes$ens_gene_id)
length(genes_use)
sum(end_drug_genes$ens_gene_id %in% early_drug_genes$ens_gene_id)
dim(rev_genes)
# obs_genes <- genes_use[!genes_use %in% rev_genes$ens_gene_id]
# length(obs_genes)



preprocess_mat <- as.matrix(counts(ts_sce[genes_use,]))
dim(preprocess_mat)

genes_use1 <- obs_genes_df %>%
  filter(gene_type=='cls_7') %>%
  pull(genes_use)
length(genes_use1)
preprocess_mat <- as.matrix(counts(ts_sce[genes_use1,]))
dim(preprocess_mat)
gene_module_df1 <- get_genes_modules(preprocess_mat=preprocess_mat, resolution=1.2)
summary(as.factor(gene_module_df1$module))
obs_genes_df1 <- data.frame(genes_use=c(gene_module_df1$id),
                           gene_type=c(paste0('cls_',gene_module_df1$module)),
                           row.names=gene_module_df1$id)
dim(obs_genes_df1)
obs_genes_df1$gene_type[1]
plttitle <- 'Genes modules 7'
p <- viz_heatmap(ts_sce, obs_genes_df1, output_dir, datatag, plttitle)
trans_genes_7 <- obs_genes_df1 %>%
  filter(gene_type %in% c('cls_1', 'cls_3')) %>%
  pull(genes_use)

repressed_genes_7 <- obs_genes_df1 %>%
  filter(!gene_type %in% c('cls_1', 'cls_3')) %>%
  pull(genes_use)

trans_genes_main <- obs_genes_df %>%
  filter(gene_type %in% c('cls_9','cls_8','cls_4','cls_12','cls_5','cls_11')) %>%
  pull(genes_use)

repressed_genes_main <- obs_genes_df %>%
  filter(gene_type %in% c('cls_3','cls_2','cls_10')) %>%
  pull(genes_use)

activated_genes_main <- obs_genes_df %>%
  filter(gene_type %in% c('cls_6','cls_1')) %>%
  pull(genes_use)

length(activated_genes_main)
length(c(repressed_genes_main, repressed_genes_7))
length(c(trans_genes_7, trans_genes_main))
write_csv_genes_modules(c(trans_genes_7, trans_genes_main), datatag, 'transient_genes', output_dir)
write_csv_genes_modules(c(repressed_genes_main, repressed_genes_7), datatag, 'repressed_genes', output_dir)
write_csv_genes_modules(activated_genes_main, datatag, 'activated_genes', output_dir)
tag <- 'transient_genes'
trans_df <- data.table::fread(paste0(output_dir, datatag,'_', tag, '.csv'))
dim(trans_df)
trans_df$gene_type <- paste0('Trans(',dim(trans_df)[1],')')

tag <- 'repressed_genes'
repressed_df <- data.table::fread(paste0(output_dir, datatag,'_', tag, '.csv'))
dim(repressed_df)
repressed_df$gene_type <- paste0('Repressed(',dim(repressed_df)[1],')')

tag <- 'activated_genes'
activated_df <- data.table::fread(paste0(output_dir, datatag,'_', tag, '.csv'))
dim(activated_df)
activated_df$gene_type <- paste0('Activated(',dim(activated_df)[1],')')
total_genes <- dplyr::bind_rows(repressed_df, activated_df, trans_df)

data.table::fwrite(total_genes, paste0(output_dir, datatag,'_total_genes_modules_act_repr_trans.csv'))
total_genes <- total_genes %>%
  dplyr::rename(genes_use=ens_gene_id)
plttitle <- 'Activated Repressed Transient genes'
p <- viz_heatmap(ts_sce, total_genes, output_dir, datatag, plttitle)
p



genes_df$genes_use <- genes_df$ens_gene_id
total_genes <- total_genes %>%
  dplyr::rename(genes_use=ens_gene_id)
plttitle <- 'Activated Repressed Transient genes'
genes_df$gene_type <- genes_df$gene_type_module
unique(genes_df$gene_type)
phm <- viz_heatmap(ts_sce, genes_df, save_dir, datatag, plttitle)
phm


preprocess_mat <- as.matrix(counts(ts_sce[trans_df$ens_gene_id,]))
dim(preprocess_mat)


gene_module_df <- get_genes_modules(preprocess_mat=preprocess_mat, resolution=0.05)
# gene_module_df <- get_genes_modules(preprocess_mat=preprocess_mat, resolution=0.1)
summary(as.factor(gene_module_df$module))
colnames(gene_module_df)
gene_module_df <- gene_module_df %>%
  # dplyr::rename(ens_gene_id=id) %>%
  dplyr::select(ens_gene_id, module)

trans_df <- trans_df %>% inner_join(gene_module_df,by=c('ens_gene_id'))
summary_trans <- trans_df %>%
  dplyr::group_by(module) %>%
  dplyr::summarise(nb_cells=n())

trans_df <- trans_df %>% inner_join(summary_trans,by=c('module'))
dim(trans_df)
trans_df$gene_type <- paste0('Trans',trans_df$module,'(',trans_df$nb_cells,')')

unique(trans_df$gene_type)



data.table::fwrite(gene_module_df, paste0(output_dir,'Rx_vs_UnRx/gene_modules_',datatag,'.csv'))
dim(gene_module_df)


obs_genes_df <- data.frame(genes_use=c(gene_module_df$id),
                           gene_type=c(paste0('cls_',gene_module_df$module)),
                           row.names=gene_module_df$id)
dim(obs_genes_df)
plttitle <- 'Activated Repressed genes'

p <- viz_heatmap(ts_sce, obs_genes_df, output_dir, datatag, plttitle)

plttitle <- 'Transient_genes'
p <- viz_heatmap(ts_sce, obs_genes_df, output_dir, datatag, plttitle)

plttitle <- 'Activated genes modules'
p <- viz_heatmap(ts_sce, obs_genes_df, output_dir, datatag, plttitle)


obs_genes_df <- data.frame(genes_use=c(gene_module_df$id, rev_genes$ens_gene_id),
                           gene_type=c(paste0('Cls',gene_module_df$module),
                                       rep('transient_genes',length(rev_genes$ens_gene_id))))
repressed_cls <- paste0('Cls',c(6,7,3,4))
activated_cls <- paste0('Cls',c(2,9,8,5,1,10))
obs_genes_df$gene_type <- ifelse(obs_genes_df$gene_type %in% repressed_cls,'Repressed',
                                 ifelse(obs_genes_df$gene_type %in% activated_cls,'Activated','Transient'))
s <- obs_genes_df %>%
  group_by(gene_type)%>%
  summarise(nb_gene=n())%>%
  mutate(gene_type_desc=paste0(gene_type,'(',nb_gene,')'))
  as.data.frame()
obs_genes_df <- obs_genes_df %>% inner_join(s, by='gene_type')%>%
  select(-gene_type)
summary(as.factor(obs_genes_df$gene_type_desc))

data.table::fwrite(obs_genes_df, paste0(output_dir,'observe_genes.csv'))
t <- gene_module_df %>%
  filter(module==11)%>%
  pull(id)
sum(t %in% rev_genes$ens_gene_id)
t
dim(obs_genes_df)

dim(rev_genes)
rev_genes$pattern[1:10]
rev_genes$gene_symbol[1:50]
summary(rev_genes$pattern)
summary(rev_genes$statSecondTest)
dim(end_drug_genes)
rev_genes_ls <- rev_genes$ens_gene_id
end_drug_genes_ls <- end_drug_genes$ens_gene_id[!end_drug_genes$ens_gene_id %in% rev_genes_ls]
early_drug_genes_ls <- early_drug_genes$ens_gene_id[!early_drug_genes$ens_gene_id %in% end_drug_genes_ls]
early_drug_genes_ls <- early_drug_genes_ls[!early_drug_genes_ls %in% rev_genes_ls]
length(intersect(end_drug_genes_ls, early_drug_genes_ls))
dim(end_drug_genes)
dim(early_drug_genes)
length(early_drug_genes_ls)
length(end_drug_genes_ls)
end_drug_genes$pattern
end_drug_genes <- end_drug_genes %>%
  dplyr::filter(!ens_gene_id %in% rev_genes_ls)%>%
  dplyr::arrange(desc(pattern))
head(end_drug_genes)




obs_gene <- end_drug_genes$ens_gene_id[2]
# obs_gene <- meta_genes %>%
#   dplyr::filter(gene_symbol==gsymb)%>%
#   dplyr::pull(genes_use)




library(pheatmap)
my_hclust_gene <- hclust(dist(tmp), method = "complete")
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 5)
# my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "cluster 1", no = "cluster 2"))

head(my_gene_col)

class(my_gene_col)
my_gene_col <- as.data.frame(my_gene_col)
my_gene_col$genes_use <- rownames(my_gene_col)
dim(my_gene_col)
dim(obs_genes_df)

obs_genes_df <- obs_genes_df %>% 
  dplyr::left_join(my_gene_col, by='genes_use')%>%
  tibble::column_to_rownames('genes_use')

p1 <- readRDS(paste0(output_dir,'clusters_genes_reversibility.rds'))  


patternRes$ens_gene_id <- rownames(patternRes)
activated_df <- activated_df %>% left_join(patternRes,by='ens_gene_id')
activated_df <- activated_df %>%
    dplyr::arrange(desc(waldStat))
head(activated_df)
gact <- 'TOP2A' #UBE2C, NUF2, TPX2, CENPF
repressed_df <- repressed_df %>% left_join(patternRes,by='ens_gene_id')
repressed_df <- repressed_df %>% 
  dplyr::arrange(desc(waldStat))
head(repressed_df[,1:3])
gpr <- 'COX6C' #PAGE2, MGST1, DCAF13, CD24
trans_df <- trans_df %>% left_join(patternRes,by='ens_gene_id')
dim(trans_df)
unique(trans_df$module)
trans_df1 <- trans_df %>%
  dplyr::filter(module==1) %>%
  dplyr::arrange(desc(waldStat),module)
head(trans_df1[,1:3])
gm1_ls <- c('ID1','SPINT2','CRABP2')

trans_df2 <- trans_df %>%
  dplyr::filter(module==2) %>%
  dplyr::arrange(desc(waldStat),module)
head(trans_df2[,1:3])
gm2 <- 'CCDC26'
gm2_ls <- c('CCDC26','NEFL','MAP1B')
trans_df3 <- trans_df %>%
  dplyr::filter(module==3) %>%
  dplyr::arrange(desc(waldStat),module)
head(trans_df3[,1:3])
gm3 <- 'HIST1H4C'
gm3_ls <- c('HIST1H4C','MYBL2','CDCA5')
# trans_df2 <- trans_df %>%
#   dplyr::filter(module==2) %>%
#   dplyr::arrange(desc(waldStat),module)
# head(trans_df2)
# gm2 <- 'CCDC26'

# plt_ls[[count]] <- p1
plt_ls <- list()
obs_genes_ls <- c('CDK1','MAEL','CDK14','FOXP1','CTSZ')
obs_genes_ls <- c('CDK14','MAEL','CTSZ')
gact_ls <- c('TOP2A', 'UBE2C', 'NUF2', 'TPX2', 'CENPF')
gpr_ls <- c('COX6C', 'PAGE2', 'MGST1', 'DCAF13', 'CD24')
obs_genes_ls <- c(gact_ls, gpr_ls, gm1_ls, gm2_ls, gm3_ls)
colnames(total_genes)
names(obs_genes_ls) <- c(rep('activated', length(gact_ls)),
                         rep('repressed', length(gpr_ls)),
                         rep('transient1', length(gm1_ls)),
                         rep('transient2', length(gm2_ls)),
                         rep('transient3', length(gm3_ls)))
figs_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/figs_v3/"
total_genes <- as.data.frame(total_genes)
# meta_genes <- total_genes %>%
#   as.data.frame() %>%
#   remove_rownames %>% 
#   dplyr::filter(gene_symbol %in% as.character(obs_genes_ls)) %>%
#   dplyr::select(genes_use, gene_symbol) %>%
#   tibble::column_to_rownames('gene_symbol')
sum(total_genes$gene_symbol %in% obs_genes_ls)

meta_genes <- total_genes[total_genes$gene_symbol %in% obs_genes_ls,]
library(tidyverse)
meta_genes <- meta_genes %>% 
  remove_rownames %>%
  dplyr::select(genes_use, gene_symbol) %>%
  tibble::column_to_rownames('gene_symbol')


viz_given_gene_exp_lineages(obs_genes_ls, meta_genes, ts_sce, figs_dir, datatag)
  

# Testing
obs_genes_ls <- c('TOP2A')
obs_genes_ls <- c('TOP2A','COX6C','ID1','CCDC26')
datatag <- 'SA609'
figs_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/figs_v3/"
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/'
total_genes <- data.table::fread(paste0(save_dir, 'SA609_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
dim(total_genes)
meta_genes <- total_genes %>%
  dplyr::filter(gene_symbol %in% obs_genes_ls)

library(tidyverse)
meta_genes <- meta_genes %>% 
  remove_rownames %>%
  dplyr::select(ens_gene_id, gene_symbol) %>%
  tibble::column_to_rownames('gene_symbol')

head(meta_genes)
avg_gene_exp <- get_average_gene_exp_per_lineage(datatag, ts_sce, obs_genes_ls, output_fn, nEstimatedPoints=100, save_data=T)
# save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/')
# ts_sce <- readRDS(paste0(save_dir,'tradeSeq_3000/', "fitGAM_out.rds"))
viz_given_gene_exp_lineages(obs_genes_ls, meta_genes, avg_gene_exp, figs_dir, datatag)

  

for(gsymb in obs_genes_ls){
  count <- count + 1
  pg <- readRDS(paste0(output_dir,gsymb,'.rds'))  
  plt_ls[[gsymb]] <- pg + theme(legend.position = 'none')
}  
pg <- readRDS(paste0(output_dir,gsymb,'.rds'))  
plg_lg <- cowplot::get_legend(pg)
plt_ls[['lg']] <- cowplot::ggdraw(plg_lg)

library(ComplexHeatmap)
p_obs_genes <- cowplot::plot_grid(plotlist = plt_ls, nrow = 1)
gb <- grid.grabExpr(draw(p1, annotation_legend_side = "bottom",
                         heatmap_legend_side = "bottom"),
                    padding = unit(c(1, 1, 2, 2), "mm"))#,merge_legend = TRUE

ptotal <- cowplot::plot_grid(gb, p_obs_genes, nrow = 2, rel_heights = c(2,1))
png(paste0(output_dir,"reverse_33genes_summary_",datatag,".png"), height = 2*650, width=2*400,res = 2*72)
print(ptotal)
dev.off()

png(paste0(output_dir,"reverse_33genes_summary_",datatag,".png"), height = 2*550, width=2*650,res = 2*72)
print(ptotal)
dev.off()

# load sce file
# Get up-regulated genes
# node_df: genes_df - genes modules only, log10 of statistical pattern test --> logFC
# extract a module of genes to test first
# edges: subset only genes in node_df
output_dir <- paste0(output_dir,'genes_network/')
dir.create(output_dir)


node_df <- as.data.frame(rev_genes)
node_df <- node_df[!duplicated(node_df$gene_symbol),]
dim(node_df)
colnames(node_df)
rownames(node_df) <- node_df$ens_gene_id
node_df$logFC <- log2(node_df$pattern)
summary(node_df$logFC)
obs_genes <- node_df$gene_symbol
# node_df$gene_symbol
rownames(sce)[1]
nex <- as.data.frame(logcounts(sce[node_df$ens_gene_id,cells_use]))
rownames(nex) <- node_df[rownames(nex),'gene_symbol']
dim(nex)
print('Load String DB...')
string_net_out_fn <- paste0(output_dir, datatag,'_string_net_stat.rds')
filtered_genes <- patternRes$gene_symbol
string_net_stat <- load_stringDB(string_net_out_fn, filtered_genes, output_dir, datatag,
             score_thrs=150, save_data=T)

string_network <- string_net_stat$string_network
dim(string_net_stat$protein_infos)
dim(string_network)
head(string_net_stat$protein_infos)
protein_infos <- string_net_stat$protein_info
rownames(protein_infos) <- protein_infos$protein_external_id
string_connection <- paste0(protein_infos[string_network$protein1,'preferred_name'],'_',
                            protein_infos[string_network$protein2,'preferred_name'])


corr_stat <- compute_correlation(nex, obs_genes, output_dir, datatag, corr_thrs=0.3)
edges_df <- get_edges_ls(corr_summary, string_connection)
construct_graph_v2(edges_df, node_df, output_dir, datatag, cluster_use='reverse_DH')




metacell <- colData(sce) %>% as.data.frame()
dim(metacell)
unique(metacell$clone)
metacell$clone <- get_unique_clone_id(metacell$clone)
unique(metacell$clone)
metacell <- metacell %>%
  dplyr::mutate(clone=replace(clone,clone=='R','A'))
unique(metacell$treatmentSt)
# metacell$treatment_desc <- get_treatment_desc(metacell$treatmentSt)
metacell$treatment_desc <- get_treatment_detail_desc(metacell$treatmentSt)

unique(metacell$treatment_desc)
table(metacell$treatment_desc,metacell$clone)
metadata <- metacell %>%
  dplyr::select(treatment_desc,clone,cell_id,treatmentSt)%>%
  dplyr::mutate(clone=paste0('clone ',clone))%>%
  dplyr::mutate(treatment_clone=paste0(treatment_desc,': clone',clone))
pseudo <- data.table::fread(paste0(save_dir,'slingshot_SA609_10_PCA_pseudotime.csv')) %>% as.data.frame()
weight <- data.table::fread(paste0(save_dir,'slingshot_SA609_10_PCA_cellWeights.csv')) %>% as.data.frame()
dim(weight)  
weight <- weight %>% inner_join(metadata,by='cell_id')
weight <- weight %>%
  dplyr::filter(Lineage2>0)
dim(weight)
cells_use <- weight$cell_id

patternRes$pattern
patternRes1 <- patternRes %>%
  filter(ens_gene_id %in% obs_genes)%>%
  dplyr::rename(logFC=pattern)
genes_df$genes_use
deg_df <- endRes12 %>%
  dplyr::select(logFC1_2, gene_symb)%>%
  dplyr::rename(logFC=logFC1_2)
deg_df <- patternRes1
dim(deg_df)
# TO DO: what is the pathway that related to RxH genes? 144 genes
gsea_out <- get_custom_pathway_results(deg_df,                      # named vector of statistical significance 
                                       desc='Up,Down regulated genes', base_name = datatag,  
                                       pathway_name=c('hallmark'), #'cosmic' or 'cisplatin_resistance', or 'metastasis'
                                       groups_use=c("Rx","UnRx"),    # vector of 2 elements, 2 group name used for DE analysis
                                       output_dir,
                                       reference_genes_set=NULL,
                                       n_top=20,
                                       ref_dif=NULL)

gsea_out$pathway
class(gsea_out)
gsea_out <- gsea_out %>%
  dplyr::filter(padj < 0.05)


save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/'
genes_df <- data.table::fread(paste0(save_dir, 'SA609_total_genes_modules_act_repr_trans.csv')) %>% as.data.frame()
dim(genes_df)
summary(as.factor(genes_df$gene_type))

meta_types <- data.frame(gene_type=unique(genes_df$gene_type),
                         gene_type_module=c('Module4','Module1','Module5','Module6',
                                            'Module2','Module3'))
genes_df <- genes_df %>% inner_join(meta_types, by='gene_type')
colnames(genes_df)
genes_df <- genes_df %>%
  dplyr::select(ens_gene_id, gene_symbol, gene_type_module, everything()) %>%
  dplyr::select(-nb_cells, -module)

head(genes_df)

data.table::fwrite(genes_df, paste0(save_dir, 'SA609_total_genes_modules_act_repr_trans_08_Dec.csv.gz'))

save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/'
genes_df <- data.table::fread(paste0(save_dir, 'SA609_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
dim(genes_df)
colnames(genes_df)


datatag <- 'SA609'
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/')
ts_sce <- readRDS(paste0(save_dir,'tradeSeq_3000/', "fitGAM_out.rds"))
plttitle <- 'Activated Repressed Transient genes'
genes_df$gene_type <- genes_df$gene_type_module
unique(genes_df$gene_type_module)

phm <- viz_heatmap(ts_sce, genes_df, save_dir, datatag, plttitle)

results_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/figs_v3/"
saveRDS(phm, paste0(results_dir,'Activated_Repressed_Transient_genes.rds'))


patternRes <- readRDS('/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/patternRes_out.rds')
dim(patternRes)
colnames(patternRes)
patternRes$ens_gene_id <- rownames(patternRes)
patternRes1 <- annotables::grch38 %>%
  # dplyr::select(gene_symbol = symbol, chr, start, end) %>%
  dplyr::select(gene_symbol = symbol, ens_gene_id = ensgene) %>%
  inner_join(patternRes, by='ens_gene_id')

# patternRes1 <- annotables::grch38 %>% 
#   # dplyr::select(gene_symbol = symbol, chr, start, end) %>% 
#   dplyr::select(gene_symbol = symbol, ens_gene_id = ensgene, entrez) %>% 
#   inner_join(patternRes, by='ens_gene_id')


deg_df <- genes_df %>% inner_join(patternRes, by='ens_gene_id')
dim(deg_df)
save_dir_pw <- paste0(save_dir,'pathway_results/')
dir.create(save_dir_pw)
datatag <- 'SA609'

gsea_out <- get_pathway_trajectory_genes(patternRes1, save_dir_pw, 
                                         desc='total_genes_modules', 
                                         base_name = datatag, pathway_name='hallmark')
gsea_out <- get_pathway_trajectory_genes(patternRes1, save_dir_pw, desc='total_genes_modules', base_name = datatag, pathway_name='kegg')
head(gsea_out)
gsea_out$pathway
gsea_out$signf_genes
deg_df1 <- deg_df %>%
  dplyr::filter(gene_type=="Repressed(168)")
gsea_out <- get_pathway_trajectory_genes(deg_df1, save_dir, desc='repressed_genes', base_name = datatag, pathway_name='hallmark')


deg_df2 <- deg_df %>%
  dplyr::filter(gene_type=="Activated(132)")
gsea_out <- get_pathway_trajectory_genes(deg_df2, save_dir, desc='activated_genes', base_name = datatag, pathway_name='hallmark')
gsea_out$pathway
gsea_out$signf_genes

gs <- c(unlist(strsplit(gsea_out$signf_genes[1],',')),
        unlist(strsplit(gsea_out$signf_genes[2],',')),
        unlist(strsplit(gsea_out$signf_genes[3],',')))
gs[[1]]
sum(gs %in% genes_df$gene_symbol)
length(gs)
deg_df2 <- deg_df %>%
  dplyr::filter(!gene_type %in% c("Activated(132)","Repressed(168)"))
dim(deg_df2)
gsea_out <- get_pathway_trajectory_genes(deg_df2, save_dir, desc='transient_genes', base_name = datatag, pathway_name='hallmark')



deg_df2$pvalue
unique(deg_df$gene_type)
dim(gsea_out)

sum(gsea_out$padj<0.05)
View(head(gsea_out))

# [1] "HALLMARK_E2F_TARGETS"     "HALLMARK_G2M_CHECKPOINT" 
# [3] "HALLMARK_MITOTIC_SPINDLE"
