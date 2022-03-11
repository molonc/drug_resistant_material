datatag <- 'SA1035'
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')
nfeatures_use <- 3000
sce <- readRDS(paste0(save_dir, datatag,'_',nfeatures_use,'_rd_sce_v2.rds'))
dim(sce)
output_dir <- paste0(save_dir, 'slingshot_output/')
plot_trajectory_clones_prevalence_SA1035(sce, output_dir, datatag)

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

