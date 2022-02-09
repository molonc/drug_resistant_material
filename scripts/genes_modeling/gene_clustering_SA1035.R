script_dir <- '/home/htran/Projects/farhia_project/rscript/genes_modeling/'
source(paste0(script_dir, "gene_utils.R"))

base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# gene_type <- 'increase'
gene_type <- 'decrease'
save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/gene_regression_',gene_type,'/')
save_dir
datatag <- 'SA1035'
obs_clones <- c('H','A_H','H_K','I','J','K')
# clone_aware <- FALSE
clone_aware <- TRUE
obs_treatment_st <- c('UT','UTT','UTTT','UTTTT','UTTTTT')

obs_clones_untreated <- c('E')
obs_untreated_st <- c('UU','UUU','UUUU','UUUUU')

# input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results'
dir.create(save_dir, recursive = T, showWarnings=F)
# Define clone R: R, clone S: H
# Read csv files of DE genes
de_x3 = read.csv(paste0(base_dir,'SA1035_rna/deg_analysis/SA1035-v6/SA1035_UTTT_H_UUUU_E/signif_genes.csv'))
de_x4 = read.csv(paste0(base_dir,'SA1035_rna/deg_analysis/SA1035-v6/SA1035_UTTTT_H_UUUUU_E/signif_genes.csv'))



# de_x3 <- get_top_up_genes(de_x3, ntop = 500)
# de_x4 <- get_top_up_genes(de_x4, ntop = 500)
# print(dim(de_x3))
sce <- readRDS(paste0(base_dir,'rnaseq_v6/SA1035-v6/SA1035_sctransform_normalized.rds')) #rnaseq_v6/SA1035-v6/total_sce_clones.rds
print(dim(sce))
print(dim(sce))
if(is.null(rowData(sce)$Symbol)){
  genes_symb_df <- read.csv(paste0(base_dir,'biodatabase/meta_genes.csv'), check.names = F, stringsAsFactors = F)
  dim(genes_symb_df)
  rownames(genes_symb_df) <- genes_symb_df$gene_ens
  rowData(sce)$Symbol <- genes_symb_df[rownames(sce),'gene_symb']
  print(rowData(sce)$Symbol[1])
}

# Select only resistant genes 
# For testing, get top genes first

if(gene_type=='increase'){
  print('option increase')
  observed_genes <- intersect(de_x3[de_x3$logFC>0,'ensembl_gene_id'], de_x4[de_x4$logFC>0,'ensembl_gene_id'])
}else{
  print('option decrease')
  observed_genes <- intersect(de_x3[de_x3$logFC<0,'ensembl_gene_id'], de_x4[de_x4$logFC<0,'ensembl_gene_id'])
}
print(summary(as.factor(de_x4[de_x4$ensembl_gene_id %in% observed_genes,'Gene_Type'])))
observed_genes <- intersect(observed_genes, rownames(sce))
# length(observed_genes)
# de_x3 <- de_x3[de_x3$ensembl_gene_id %in% observed_genes,]
# de_x3 <- de_x3 %>%
#   dplyr::filter(ensembl_gene_id %in% observed_genes)
# print(dim(de_x3))
# de_x3 <- get_top_up_genes(de_x3, ntop = 500)
# observed_genes <- as.character(de_x3$ensembl_gene_id)

print(paste0("Number of observed genes: ", length(observed_genes)))
# observed_sce <- sce[observed_genes,]
# rownames(sce)[1:3]
# unique(observed_sce$clone)
# dim(observed_sce)
# summary(as.factor(observed_sce$treatmentSt))



meta_cells_df <- as.data.frame(colData(sce))
# meta_cells_df$cell_id <- rownames(meta_cells_df)
# meta_cells_df$treatmentSt <- get_treatment_status(meta_cells_df$series)
meta_cells_df$treatmentSt[1:3]
meta_cells_df <- meta_cells_df %>%
  dplyr::select(cell_id, treatmentSt, timepoint, clone) %>%
  dplyr::mutate(treatment_status=paste0(timepoint,'_',treatmentSt))
# colnames(meta_cells_df)
# View(head(meta_cells_df[,c("series","timepoint","clone")]))

if(clone_aware){
  print("With clone aware")
  meta_cells_untreated <- meta_cells_df %>%
    dplyr::filter(clone %in% obs_clones_untreated & treatmentSt %in% obs_untreated_st)
  print(dim(meta_cells_untreated))
  meta_cells_df <- meta_cells_df %>%
    dplyr::filter(clone %in% obs_clones & treatmentSt %in% obs_treatment_st)
}else{
  print("Without clone aware")
  meta_cells_untreated <- meta_cells_df %>%
    dplyr::filter(treatmentSt %in% obs_untreated_st)
  print(dim(meta_cells_untreated))
  meta_cells_df <- meta_cells_df %>%
    dplyr::filter(treatmentSt %in% obs_treatment_st)
}
print(dim(meta_cells_df))
table(meta_cells_df$clone, meta_cells_df$treatmentSt)
print(dim(meta_cells_untreated))
table(meta_cells_untreated$clone, meta_cells_untreated$treatmentSt)
# meta_cells_df <- downsample_cells(meta_cells_df, downsample_ratio=0.4,
#                              thres_small_clone=500)
# colnames(meta_cells_df)

# meta_cells_df$timepoint[1:3]



meta_cells_df <- meta_cells_df %>%
                 dplyr::filter(treatmentSt!='UT')  # only 2 cells at this time points, remove from list

# observed_sce <- observed_sce[,meta_cells_df$cell_id]
observed_sce_untreated <- sce[observed_genes, meta_cells_untreated$cell_id]
observed_sce <- sce[observed_genes, meta_cells_df$cell_id]

# View(rowData(observed_sce))
# norm_data <- logcounts(observed_sce)
# norm_data <- as.data.frame(as.matrix(norm_data))
# norm_data$gene <- rowData(observed_sce)$Symbol
# dim(norm_data)

# meta_genes <- data.frame(gene=rowData(observed_sce)$Symbol, ens_gene=rownames(observed_sce), stringsAsFactors=F)
# rownames(de_x3) <- de_x3$ensembl_gene_id
# de_x3 <- de_x3[meta_genes$ens_gene,c("logFC","ensembl_gene_id")]
# dim(de_x3)
# rownames(de_x4) <- de_x4$ensembl_gene_id
# de_x4 <- de_x4[meta_genes$ens_gene,c("logFC","ensembl_gene_id")]
# dim(de_x4)
# colnames(de_x3)[which(names(de_x3)=='logFC')] <- 'logFC_UTTT_H_UUUU_E'
# colnames(de_x4)[which(names(de_x4)=='logFC')] <- 'logFC_UTTTT_H_UUUUU_E'


meta_genes <- data.frame(gene=rowData(observed_sce)$Symbol,ens_gene=rownames(observed_sce), stringsAsFactors=F)
# meta_genes <- meta_genes %>% left_join(genes_symb_df, by=c('ens_gene'='gene_ens'))
de_x3 <- de_x3 %>%
  dplyr::filter(ensembl_gene_id %in% meta_genes$ens_gene) %>%
  dplyr::select(logFC, ensembl_gene_id) %>%
  dplyr::rename(logFC_UTTT_H_UUUU_E=logFC)

de_x4 <- de_x4 %>%
  dplyr::filter(ensembl_gene_id %in% meta_genes$ens_gene) %>%
  dplyr::select(logFC, ensembl_gene_id, Gene_Type) %>%
  dplyr::rename(logFC_UTTTT_H_UUUUU_E=logFC)
meta_genes <- meta_genes %>% left_join(de_x3, by=c('ens_gene'='ensembl_gene_id'))
meta_genes <- meta_genes %>% left_join(de_x4, by=c('ens_gene'='ensembl_gene_id'))
colnames(meta_genes)

# Input data to WGCNA
tag <- ''
norm_data <- logcounts(observed_sce)
datExpr <- as.data.frame(as.matrix(t(norm_data))) # # now cells are rows and genes are columns
res <- remove_outliers(datExpr, save_dir, datatag, tag)
datExpr <- res$datExpr
print(dim(datExpr))
norm_data <- as.data.frame(as.matrix(norm_data))
norm_data[norm_data==0] <- NA
norm_data$gene <- rowData(observed_sce)$Symbol
dim(norm_data)

norm_data_untreated <- logcounts(observed_sce_untreated)
norm_data_untreated <- as.data.frame(as.matrix(norm_data_untreated))
norm_data_untreated[norm_data_untreated==0] <- NA
norm_data_untreated$gene <- rowData(observed_sce_untreated)$Symbol
dim(norm_data_untreated)



metacells_df <- encode_metadata(meta_cells_df, save_dir, datatag)
table(rownames(metacells_df)==rownames(datExpr))

select_soft_threshold_power(datExpr, save_dir, option=1)
softPowerThres <- 10 #10 increase, 10 decrease
minClusterSz <- 20
results_mtx <- compute_adjacency_correlation_mtx(datExpr, save_dir, softPowerThres, 1, T)

results <- generate_cluster_gene_tree_v1(meta_genes, datExpr, results_mtx$diss_mtx, save_dir, softPowerThres, minClusterSz, T)  
# results <- generate_cluster_gene_tree_v2(datExpr, results_mtx$diss_mtx, save_dir, minClusterSz)  
# correlate_genes_with_treatment_effect(results, datExpr, metacells_df, save_dir)

# Increase
# gene_cluster <- results$gene_cluster
# gene_cluster[gene_cluster$cluster %in% c('turquoise','yellow'),'cluster'] <- 'yellow'
# gene_cluster[gene_cluster$cluster %in% c('brown'),'cluster'] <- 'green'

# Decrease
gene_cluster <- results$gene_cluster
gene_cluster[gene_cluster$cluster %in% c('grey','brown'),'cluster'] <- 'blue'
gene_cluster[gene_cluster$cluster %in% c('yellow'),'cluster'] <- 'green'

# summary(as.factor(gene_cluster$cluster))
write.csv(gene_cluster, paste0(save_dir,'genes_clusters.csv'), quote=F, row.names = F)
results$gene_cluster <- gene_cluster
save(results, file= paste0(save_dir, "Network_v1.RData"))

if(gene_type=='increase'){
  plttitle <- paste0(datatag,': up-regulated DEG in Rx vs UnRx')
}else{
  plttitle <- paste0(datatag,': down-regulated DEG in Rx vs UnRx')
}
plttitle
#results$gene_cluster
plot_genes_exp(gene_cluster, meta_genes, norm_data, norm_data_untreated,meta_cells_untreated,
               meta_cells_df, obs_clones, save_dir, datatag, clone_aware=TRUE, plttitle)

# plot_genes_exp(results$gene_cluster, meta_genes, norm_data, norm_data_untreated,meta_cells_untreated,
#                meta_cells_df, obs_clones, save_dir, datatag, clone_aware=TRUE)
# plot_actual_genes_exp(results$gene_cluster, meta_genes, norm_data, norm_data_untreated,meta_cells_untreated,
#                       meta_cells_df, obs_clones, save_dir, datatag, clone_aware=TRUE)


clusters_ls <- unique(results$gene_cluster$cluster)
for(cluster_use in clusters_ls){
  net <- extract_genes_correlation_network(results_mtx, results, datExpr, save_dir, datatag, softPowerThres, cluster_use)
  construct_graph(net$edges_df, net$node_df, save_dir, datatag, cluster_use)
}

load(paste0(save_dir, "Network_v1.RData"))  # load results objects
tag <- ''
datExpr <- data.table::fread(paste0(save_dir, datatag,'_',tag,'_datExpr.csv.gz'))
datExpr <- as.data.frame(datExpr)
rownames(datExpr) <- datExpr$V1
datExpr$V1 <- NULL
metacells_df <- read.csv(paste0(save_dir, datatag,'_meta_cells.csv'), check.names = F, stringsAsFactors = F)

rownames(metacells_df) <- metacells_df$cell_id
metacells_df$cell_id <- NULL
correlate_genes_with_treatment_effect(results, datExpr, metacells_df, save_dir, module = "turquoise")



# cluster_use <- 'yellow'
# net <- extract_genes_correlation_network(results_mtx, results, datExpr, save_dir, datatag, softPowerThres, cluster_use)
# construct_graph(net$edges_df, net$node_df, save_dir, datatag, cluster_use)
# 
# cluster_use <- 'black'
# net <- extract_genes_correlation_network(results_mtx, results, datExpr, save_dir, datatag, softPowerThres, cluster_use)
# construct_graph(net$edges_df, net$node_df, save_dir, datatag, cluster_use)
# 
# cluster_use <- 'green'
# net <- extract_genes_correlation_network(results_mtx, results, datExpr, save_dir, datatag, softPowerThres, cluster_use)
# construct_graph(net$edges_df, net$node_df, save_dir, datatag, cluster_use)

# gene_clustering(meta_genes, norm_data, meta_cells_df, obs_clones, 
#                 save_dir, datatag, TRUE)


