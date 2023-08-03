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
  library(clusterExperiment)
})



datatag <- 'SA1035' # Patient 6
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')
output_dir <- paste0(save_dir,'tradeseq_v2/')

if(!file.exists(output_dir)) dir.create(output_dir)
nfeatures_use <- 3000
norm_sce <- readRDS(paste0(save_dir, datatag,'_',nfeatures_use,'_rd_sce_v2.rds'))
dim(norm_sce)

var_genes_df <- data.table::fread(paste0(save_dir, "SA1035_3000_hvg_genes.csv")) %>% as.data.frame()
# View(head(var_genes_df))
rownames(norm_sce)[1]
norm_sce <- norm_sce[var_genes_df$bc[1:1000],] # get results first
counts <- as.matrix(counts(norm_sce))
dim(counts)

crv <- readRDS(paste0(save_dir, "slingshot_output/slingshot_pseudotime_SA1035_0_PCA_crv.rds"))
set.seed(7)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)

set.seed(5)
# icMat <- evaluateK(counts = counts, sds = crv, k = 3:10, 
#                    nGenes = 200, verbose = T)

run_tradeSeq(counts, pseudotime, cellWeights, output_dir)

load_tradeSeq_results <- function(){
  ts_sce <- readRDS(paste0(output_dir, "fitGAM_out.rds"))
  startRes <- readRDS(paste0(output_dir, "startRes_out.rds"))
  endRes <- readRDS(paste0(output_dir, "endRes_out.rds"))
  patternRes <- readRDS(paste0(output_dir, "patternRes_out.rds"))
  earlyDERes <- readRDS(paste0(output_dir, "earlyDERes_out.rds"))
}

datatag <- 'SA1035'
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')
nfeatures_use <- 3000
sce <- readRDS(paste0(save_dir, datatag,'_',nfeatures_use,'_rd_sce_v2.rds'))
dim(sce)
output_dir <- paste0(save_dir, 'slingshot_output/')
plot_trajectory_clones_prevalence_SA1035(sce, output_dir, datatag)

output_stat_dir <- paste0(output_dir,'results/')
dir.create(output_stat_dir)
end_df <- get_square_rank(patternRes, endRes, output_stat_dir, base_name='patternTest_diffEndTest')
early_df <- get_square_rank(patternRes, earlyDERes, output_stat_dir, base_name='patternTest_earlyDETest')
dim(end_df)
dim(early_df)
thrs_stat <- 200
statSecondTest_thrs <- 40
summary(end_df$statSecondTest)
end_df <- end_df %>%
  filter(pattern>thrs_stat & statSecondTest>statSecondTest_thrs)%>%
  dplyr::arrange(desc(pattern))
early_df <- early_df %>%
  filter(pattern>thrs_stat & statSecondTest>statSecondTest_thrs)%>%
  dplyr::arrange(desc(pattern))
dim(end_df)
dim(early_df)

dim(rev_genes)
genes_use <- unique(c(end_df$ens_gene_id, early_df$ens_gene_id))
length(genes_use)
                            

preprocess_mat <- as.matrix(counts(ts_sce[genes_use,]))
dim(preprocess_mat)

gene_module_df <- get_genes_modules(preprocess_mat=preprocess_mat, resolution=0.1)
summary(as.factor(gene_module_df$module))
data.table::fwrite(gene_module_df, paste0(output_stat_dir,'gene_modules_',datatag,'_total.csv'))
# gene_module_df$id[1]
obs_genes_df <- data.frame(genes_use=gene_module_df$id,
                           gene_type=paste0('cls_',gene_module_df$module), 
                           row.names=gene_module_df$id)
sig_firstts <- gene_module_df %>%
  filter(module %in% c(8,7))
dim(sig_firstts)
data.table::fwrite(sig_firstts, paste0(output_stat_dir,'gene_modules_',datatag,'_first_treatment.csv'))

transient_df <- gene_module_df %>%
  filter(!module %in% c(8,7))
dim(transient_df)
data.table::fwrite(transient_df, paste0(output_stat_dir,'gene_modules_',datatag,'_transient_genes.csv'))



plttitle <- 'total_genes'
p <- viz_heatmap(ts_sce, obs_genes_df, output_stat_dir, datatag, plttitle)

obs_modules <- c(8,7)
plttitle <- 'First treatment genes'
obs_genes_df1 <- obs_genes_df %>%
  filter(gene_type %in% paste0('cls_',obs_modules))
dim(obs_genes_df1)

p <- viz_heatmap(ts_sce, obs_genes_df1, 
                 output_stat_dir, datatag, plttitle)

plttitle <- 'Transient genes'
obs_genes_df2 <- obs_genes_df %>%
  filter(!gene_type %in% paste0('cls_',obs_modules))
dim(obs_genes_df2)

p <- viz_heatmap(ts_sce, obs_genes_df2, 
                 output_stat_dir, datatag, plttitle)

datatag <- 'SA1035'
output_stat_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/tradeseq_v2/results/'
transient_df <- data.table::fread(paste0(output_stat_dir,'gene_modules_',datatag,'_transient_genes.csv'))
dim(transient_df)
head(transient_df)
total_df <- data.table::fread(paste0(output_stat_dir,'gene_modules_',datatag,'_total.csv'))
head(total_df)

gap <- readRDS(paste0(results_dir,'total_signif_genes_hm.rds'))
gap
meta_genes <- gap@matrix_param$row_split
class(meta_genes)
head(meta_genes)
meta_genes <- meta_genes %>%
  select(gene_type)
meta_genes$ens_gene_id <- rownames(meta_genes)
data.table::fwrite(meta_genes, paste0(output_stat_dir, 'genes_clusters_total_SA1035.csv'))
output_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/tradeseq_v2/"
ts_sce <- readRDS(paste0(output_dir, "fitGAM_out.rds"))



total_genes <- meta_genes %>% left_join(annots, by=c('ens_gene_id'))
data.table::fwrite(total_genes, paste0(output_stat_dir, datatag,'_total_genes_modules_act_repr_trans.csv'))
total_genes <- total_genes %>%
  dplyr::rename(genes_use=ens_gene_id)
plttitle <- 'Activated Repressed Transient genes'
total_genes$gene_type[1]
total_genes$gene_type <- gsub('Transient','Trans',total_genes$gene_type)
p <- viz_heatmap(ts_sce, total_genes, save_figs_dir, datatag, plttitle)
p

unique(total_genes$gene_type)
total_genes$genes_use[1]
patternRes$ens_gene_id <- rownames(patternRes)
total_genes$ens_gene_id <- total_genes$genes_use
firstts_genes <- total_genes %>% 
  dplyr::filter(grepl('1Rx',gene_type)) %>% 
  left_join(patternRes,by='ens_gene_id')%>%
  dplyr::arrange(desc(waldStat))%>%
  dplyr::pull(ens_gene_id)
gfirst_ls <- firstts_genes[1:3]

trans_df <- total_genes %>% 
  left_join(patternRes,by='ens_gene_id') 

dim(trans_df)
unique(trans_df$gene_type)
trans_genes1 <- trans_df %>%
  dplyr::filter(grepl('Trans1',gene_type)) %>%
  dplyr::arrange(desc(waldStat))%>%
  dplyr::pull(ens_gene_id)

gm1_ls <- trans_genes1[1:2]

trans_genes2 <- trans_df %>%
  dplyr::filter(grepl('Trans2',gene_type)) %>%
  dplyr::arrange(desc(waldStat))%>%
  dplyr::pull(ens_gene_id)

gm2_ls <- trans_genes2[1:2]

trans_genes3 <- trans_df %>%
  dplyr::filter(grepl('Trans3',gene_type)) %>%
  dplyr::arrange(desc(waldStat))%>%
  dplyr::pull(ens_gene_id)

gm3_ls <- trans_genes3[1:2]

trans_genes4 <- trans_df %>%
  dplyr::filter(grepl('Trans4',gene_type)) %>%
  dplyr::arrange(desc(waldStat))%>%
  dplyr::pull(ens_gene_id)

gm4_ls <- trans_genes4[1:2]

# plt_ls[[count]] <- p1
plt_ls <- list()

obs_genes_ls <- c(gfirst_ls, gm1_ls, gm2_ls, gm3_ls, gm4_ls)
colnames(total_genes)
names(obs_genes_ls) <- c(rep('first treatment', length(gfirst_ls)),
                         rep('transient1', length(gm1_ls)),
                         rep('transient2', length(gm2_ls)),
                         rep('transient3', length(gm3_ls)),
                         rep('transient4', length(gm3_ls)))



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
viz_given_gene_exp_lineages(obs_genes_ls1, meta_genes, ts_sce, save_figs_dir, datatag)

# preprocess_mat <- as.matrix(counts(ts_sce[transient_df$id,]))
# dim(preprocess_mat)
# gene_module_df <- get_genes_modules(preprocess_mat=preprocess_mat, resolution=0.03)
# summary(as.factor(gene_module_df$module))
# class(gene_module_df)
# 
# gene_module_df <- gene_module_df %>%
#   as.data.frame() %>%
#   dplyr::rename(gene_cluster=module)%>%
#   select(id, gene_cluster)


transient_df <- transient_df %>% left_join(gene_module_df,by='id')
unique(transient_df$gene_cluster)
transient_df$gene_cluster <- paste0('Transient',transient_df$gene_cluster)
sig_firstts <- data.table::fread(paste0(output_stat_dir,'gene_modules_',datatag,'_first_treatment.csv'))
dim(sig_firstts)
sig_firstts$gene_cluster <- '1Rx'

genes_df <- dplyr::bind_rows(transient_df, sig_firstts)
dim(genes_df)
obs_genes_df <- genes_df %>%
  dplyr::rename(genes_use=id, gene_type=gene_cluster)%>%
  select(genes_use, gene_type)
rownames(obs_genes_df) <- obs_genes_df$genes_use
unique(obs_genes_df$gene_type)
# obs_genes_df$gene_type <- factor(obs_genes_df$gene_type, levels=c("1Rx","Transient_1",
#                                  "Transient_2","Transient_3","Transient_4"))
obs_genes_stat <- obs_genes_df %>%
  dplyr::group_by(gene_type) %>%
  dplyr::summarise(nb_genes=n())
obs_genes_df <- obs_genes_df %>% inner_join(obs_genes_stat,by='gene_type')
obs_genes_df$gene_type <- paste0(obs_genes_df$gene_type,'(',obs_genes_df$nb_genes,')')
unique(obs_genes_df$gene_type)
obs_genes_df$gene_type <- gsub('Transient_','Transient',obs_genes_df$gene_type)
plttitle <- 'total_signif_genes_hm'
p <- viz_heatmap(ts_sce, obs_genes_df, output_stat_dir, datatag, plttitle)



# Change the name of the modules
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/tradeseq_v2/results/'
genes_df <- data.table::fread(paste0(save_dir, 'SA1035_total_genes_modules_act_repr_trans.csv')) %>% as.data.frame()
dim(genes_df)
summary(as.factor(genes_df$gene_type))

meta_types <- data.frame(gene_type=unique(genes_df$gene_type),
                         gene_type_module=c('Module1','Module3','Module4',
                                            'Module2','Module5'))
genes_df <- genes_df %>% inner_join(meta_types, by='gene_type')
colnames(genes_df)
genes_df <- genes_df %>%
  dplyr::select(ens_gene_id, gene_symbol, gene_type_module, everything())


head(genes_df)

data.table::fwrite(genes_df, paste0(save_dir, 'SA1035_total_genes_modules_act_repr_trans_08_Dec.csv.gz'))

genes_df <- data.table::fread(paste0(save_dir, 'SA1035_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
dim(genes_df)
patternRes <- readRDS('/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/tradeseq_v2/patternRes_out.rds')
dim(patternRes)
patternRes$ens_gene_id <- rownames(patternRes)
save_dir_pw <- paste0(save_dir,'pathway_results/')
dir.create(save_dir_pw)
datatag <- 'SA1035'
patternRes1 <- annotables::grch38 %>%
  # dplyr::select(gene_symbol = symbol, chr, start, end) %>%
  dplyr::select(gene_symbol = symbol, ens_gene_id = ensgene) %>%
  inner_join(patternRes, by='ens_gene_id')

deg_df <- genes_df %>% inner_join(patternRes, by='ens_gene_id')
dim(deg_df)
gsea_out <- get_pathway_trajectory_genes(patternRes1, save_dir_pw, desc='total_genes_modules', base_name = datatag, pathway_name='hallmark')
gsea_out <- get_pathway_trajectory_genes(patternRes1, save_dir_pw, desc='total_genes_modules', base_name = datatag, pathway_name='kegg')
dim(gsea_out)
gsea_out$pathway
paste0(save_dir, pathway_name, "_signf_pathways_",base_name,"_",desc,".csv")
gsea_out <- data.table::fread(paste0(save_dir, 'pathway_results/')) %>% as.data.frame()



## Visualize specific genes expr
datatag <- 'SA1035'
figs_dir <- paste0("/home/htran/storage/datasets/drug_resistance/rna_results/",datatag,"_rna/slingshot_trajectory/figs_v3/")
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
total_genes <- data.table::fread(paste0(save_dir, datatag,'_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
dim(total_genes)
obs_genes_ls <- total_genes$ens_gene_id
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/'
ts_sce <- readRDS(paste0(input_dir,'tradeseq_v2/', "fitGAM_out.rds"))
dim(ts_sce)
output_fn <- paste0(input_dir,'tradeseq_v2/','sigf_gene_exp.csv.gz')
avg_gene_exp <- get_average_gene_exp_per_lineage(datatag, ts_sce, obs_genes_ls, output_fn, nEstimatedPoints=100, save_data=T)

## Loading from file
avg_gene_exp <- data.table::fread(output_fn) %>% as.data.frame()
# save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/')
dim(avg_gene_exp)
obs_genes_ls <- c('S100A6','CST3','NFKBIA','CSTB')
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
datatag <- 'SA1035'
genes_df <- data.table::fread(paste0(save_dir, datatag,'_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
dim(genes_df)
colnames(genes_df)
# summary(as.factor(genes_df$gene_type_module))
genes_df <- genes_df %>%
  dplyr::select(-gene_type) %>%
  dplyr::rename(gene_type=gene_type_module)
summary(as.factor(genes_df$gene_type))

plttitle <- 'Activated Repressed Transient genes'
save_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/"
# ts_sce <- readRDS(paste0(save_dir,'tradeSeq_v2/', "fitGAM_out.rds"))
dim(ts_sce)
obs_lineages <- c(1,2,3,4) ## lineages list that you want to display
phm <- viz_heatmap(ts_sce, genes_df, save_dir, datatag, plttitle, obs_lineages)
output_dir <- save_dir
exp_mtx <- data.table::fread(paste0(output_dir,"mtx_hm.csv.gz")) %>% as.data.frame()
obs_genes_df <- data.table::fread(paste0(output_dir,"obs_genes_hm.csv.gz")) %>% as.data.frame()
obs_cells_df <- data.table::fread(paste0(output_dir,"obs_cells_hm.csv.gz")) %>% as.data.frame()
rownames(exp_mtx) <- exp_mtx$ens_gene_id
exp_mtx$ens_gene_id <- NULL
dim(obs_cells_df)
head(obs_cells_df)
dim(exp_mtx)
dim(obs_genes_df)
dim(obs_cells_df)
input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/cis_trans_lineages/')
cistrans_anno <- data.table::fread(paste0(input_dir, datatag,'_genes_cis_trans_lineages.csv')) %>% as.data.frame()
dim(cistrans_anno)
# head(cistrans_anno)
rownames(cistrans_anno) <- NULL
cistrans_anno$gene_type <- paste0("in ",cistrans_anno$gene_type)
summary(as.factor(cistrans_anno$gene_type)). ## To Do: check input data here
meta_clone_lg <- data.table::fread(paste0(input_dir,datatag,'_meta_clone_lineages.csv')) %>% as.data.frame()
meta_clone_lg <- meta_clone_lg[!duplicated(meta_clone_lg$lineage),]
meta_clone_lg
meta_clone_lg$lineage_desc <- gsub(' ','',meta_clone_lg$lineage_desc)
meta_clone_lg <- meta_clone_lg %>%  # less important, excluded from analysis here
  dplyr::filter(lineage!='Lineage5')
# sum(obs_genes_df$ens_gene_id==rownames(exp_mtx))
unique(obs_genes_df$gene_type)

gap <- viz_genes_exp_lineages_cistrans_anno_hm(as.matrix(exp_mtx), obs_genes_df, obs_cells_df,
                                               cistrans_anno, meta_clone_lg,
                                               output_dir, plttitle)
gap
