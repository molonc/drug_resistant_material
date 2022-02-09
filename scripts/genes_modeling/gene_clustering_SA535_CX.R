script_dir <- '/home/htran/Projects/farhia_project/rscript/genes_modeling/'
source(paste0(script_dir, "gene_utils.R"))

base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
gene_type <- 'increase'
# gene_type <- 'decrease'
save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/gene_regression_',gene_type,'_cx/')
save_dir

datatag <- 'SA535'
obs_clones <- c('R','U')
# clone_aware <- FALSE
clone_aware <- TRUE
obs_treatment_st <- c('UX','UXX','UXXX','UXXXX','UXXXXX')
# input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results'
dir.create(save_dir, recursive = T, showWarnings = F)
# Define clone R: R, clone S: H

obs_clones_untreated <- c('J')
obs_untreated_st <- c('UU','UUU','UUUU','UUUUU')

# Read csv files of DE genes
de_x3 = read.csv(paste0(base_dir,'SA535_total_rna_v2/SA535-v6/SA535_UXXXX_U_UUUUU_J/signif_genes.csv'))
# de_x4 = read.csv(paste0(base_dir,'SA535_total_rna_v2/SA535-v6/SA535_UUTTT_T_UUUUU_J/signif_genes.csv'))
dim(de_x3)
# dim(de_x4)
de_x3 <- de_x3 %>%
       dplyr::filter(!(abs(logFC)<0.5 & grepl('In_trans',Gene_Type)))

# rownames(de_x3) <- de_x3$ensembl_gene_id
# rv_genes <- sapply(de_x3$ensembl_gene_id, function(g) 
#   {
#    if(abs(de_x3[g,'logFC'])<0.5 & grepl('In_trans',de_x3[g,'Gene_Type'])){
#      return(TRUE)
#    }else{
#      return(FALSE)
#    }
#   })
# length(rv_genes)
# sum(rv_genes)
# de_x3 <- de_x3[!rv_genes,]
# de_x3 <- get_top_up_genes(de_x3, ntop = 500)
# de_x4 <- get_top_up_genes(de_x4, ntop = 500)
# print(dim(de_x3))
# total_clones <- read.csv(paste0(base_dir,'rnaseq_v6/SA535-v6/total_clones.csv'))
# table(total_clones$clone, total_clones$treatmentSt)

sce <- readRDS(paste0(base_dir,'rnaseq_v6/SA535-v6/SA535_CX5461_sctransform_normalized.rds'))  #sce_cx_clones.rds
# rownames(sce) <- rowData(sce)$ensgene
print(dim(sce))
if(is.null(rowData(sce)$Symbol)){
  genes_symb_df <- read.csv(paste0(base_dir,'biodatabase/meta_genes.csv'), check.names = F, stringsAsFactors = F)
  dim(genes_symb_df)
  rownames(genes_symb_df) <- genes_symb_df$gene_ens
  rowData(sce)$Symbol <- genes_symb_df[rownames(sce),'gene_symb']
}


# sce$clone <- ifelse(sce$clone=='TRUE','T',sce$clone)
# sum(sce$clone=='TRUE')
print(summary(as.factor(sce$clone)))
# Select only resistant genes 
# For testing, get top genes first
# observed_genes <- intersect(de_x3$ensembl_gene_id, de_x4$ensembl_gene_id)
# observed_genes <- de_x3$ensembl_gene_id
if(gene_type=='increase'){
  print('option increase')
  # observed_genes <- intersect(de_x3[de_x3$logFC>0,'ensembl_gene_id'], de_x4[de_x4$logFC>0,'ensembl_gene_id'])
  observed_genes <- de_x3[de_x3$logFC>0,'ensembl_gene_id']
}else{
  print('option decrease')
  observed_genes <- de_x3[de_x3$logFC<0,'ensembl_gene_id']
  # observed_genes <- intersect(de_x3[de_x3$logFC<0,'ensembl_gene_id'], de_x4[de_x4$logFC<0,'ensembl_gene_id'])
}

observed_genes <- intersect(observed_genes, rownames(sce))
length(observed_genes)
# de_x3 <- de_x3[de_x3$ensembl_gene_id %in% observed_genes,]
# de_x3 <- de_x3 %>%
#   dplyr::filter(ensembl_gene_id %in% observed_genes)
# 
# logFC_max_thrs <- 3.5
# de_x3 <- de_x3 %>%
#   dplyr::filter(logFC <= logFC_max_thrs)  # remove outlier genes that scale down analysis
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
print(unique(meta_cells_df$treatmentSt))
print(summary(as.factor(meta_cells_df$clone)))
meta_cells_df <- meta_cells_df %>%
  dplyr::select(cell_id, treatmentSt, timepoint, clone) %>%
  dplyr::mutate(treatment_status=paste0(timepoint,'_',treatmentSt))

if(clone_aware){
  print("With clone aware")
  meta_cells_untreated <- meta_cells_df %>%
    dplyr::filter(clone %in% obs_clones_untreated & treatmentSt %in% obs_untreated_st)#  
  
  
  meta_cells_df <- meta_cells_df %>%
    dplyr::filter(clone %in% obs_clones & treatmentSt %in% obs_treatment_st)
  
}else{
  print("Without clone aware")
  meta_cells_untreated <- meta_cells_df %>%
    dplyr::filter(treatmentSt %in% obs_untreated_st)
  
  meta_cells_df <- meta_cells_df %>%
    dplyr::filter(treatmentSt %in% obs_treatment_st)
}
print(summary(meta_cells_untreated$treatmentSt))
print(summary(meta_cells_df$treatmentSt))
print(dim(meta_cells_df))
print(dim(meta_cells_untreated))
table(meta_cells_df$clone, meta_cells_df$treatmentSt)
table(meta_cells_untreated$clone, meta_cells_untreated$treatmentSt)


observed_sce_untreated <- sce[observed_genes, meta_cells_untreated$cell_id]
observed_sce <- sce[observed_genes, meta_cells_df$cell_id]
dim(observed_sce)
dim(observed_sce_untreated)
# View(rowData(observed_sce))
dim(observed_sce)
# norm_data <- logcounts(observed_sce)
# norm_data <- as.data.frame(as.matrix(norm_data))
# norm_data$gene <- rowData(observed_sce)$Symbol
# dim(norm_data)


meta_genes <- data.frame(gene=rowData(observed_sce)$Symbol, ens_gene=rownames(observed_sce), stringsAsFactors=F)
rownames(de_x3) <- de_x3$ensembl_gene_id
de_x3 <- de_x3[meta_genes$ens_gene,c("logFC","ensembl_gene_id","Gene_Type")]
dim(de_x3)
colnames(de_x3)[which(names(de_x3)=='logFC')] <- 'logFC_UXXXX_U_UUUUU_J'


meta_genes <- meta_genes %>% left_join(de_x3, by=c('ens_gene'='ensembl_gene_id'))
dim(meta_genes)
colnames(meta_genes)
summary(as.factor(meta_genes$Gene_Type))
# gene_clustering(meta_genes, norm_data, meta_cells_df, obs_clones, 
#                 save_dir, datatag, TRUE)


# Input data to WGCNA
norm_data <- logcounts(observed_sce)
datExpr <- as.data.frame(as.matrix(t(norm_data))) # # now cells are rows and genes are columns
tag <- 'CX5461'
res <- remove_outliers(datExpr, save_dir, datatag, tag)
datExpr <- res$datExpr
print(dim(datExpr))
norm_data <- as.data.frame(as.matrix(norm_data))
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
softPowerThres <- 9 #10 for increase, 9 for decrease
minClusterSz <- 20
results_mtx <- compute_adjacency_correlation_mtx(datExpr, save_dir, softPowerThres, 1, T)

results <- generate_cluster_gene_tree_v1(meta_genes, datExpr, results_mtx$diss_mtx, save_dir, softPowerThres, minClusterSz, T)  
# results <- generate_cluster_gene_tree_v2(datExpr, results_mtx$diss_mtx, save_dir, minClusterSz)
# correlate_genes_with_treatment_effect(results, datExpr, metacells_df, save_dir)

gene_cluster <- results$gene_cluster
gene_cluster[gene_cluster$cluster %in% c('grey','green'),'cluster'] <- 'blue'
gene_cluster[gene_cluster$cluster %in% c('turquoise'),'cluster'] <- 'brown'
summary(as.factor(gene_cluster$cluster))
write.csv(gene_cluster, paste0(save_dir,'genes_clusters.csv'), quote=F, row.names = F)
results$gene_cluster <- gene_cluster
save(results, file= paste0(save_dir, "Network_v1.RData"))

results$gene_cluster
if(gene_type=='increase'){
  plttitle <- paste0(datatag,' CX5461: up-regulated DEG in Rx vs UnRx')
}else{
  plttitle <- paste0(datatag,' CX5461: down-regulated DEG in Rx vs UnRx')
}
plttitle
plot_genes_exp(gene_cluster, meta_genes, norm_data, norm_data_untreated, meta_cells_untreated,
               meta_cells_df, obs_clones, save_dir, datatag, clone_aware=TRUE, plttitle)



# rm(results)
# rm(datExpr)
# rm(metacells_df)
load(paste0(save_dir, "Network_v1.RData"))  # load results objects
tag <- 'CX5461'
datExpr <- data.table::fread(paste0(save_dir, datatag,'_',tag,'_datExpr.csv.gz'))
datExpr <- as.data.frame(datExpr)
rownames(datExpr) <- datExpr$V1
datExpr$V1 <- NULL
metacells_df <- read.csv(paste0(save_dir, datatag,'_meta_cells.csv'), check.names = F, stringsAsFactors = F)

rownames(metacells_df) <- metacells_df$cell_id
metacells_df$cell_id <- NULL
correlate_genes_with_treatment_effect(results, datExpr, metacells_df, save_dir)
