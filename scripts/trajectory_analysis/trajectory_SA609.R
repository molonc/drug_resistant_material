# TO DO
# Get cells in treated, clone R 
# Plot trajectory 
# Plot DE genes with time as conditions 
# Plot autocorrelation genes 

# Get cells in treated, and untreated clone R, H
# Plot trajectory 
# Plot DE genes with clones and late time points as conditions 
# Plot autocorrelation genes, optional 
# Get genes modules 
suppressPackageStartupMessages({
  library(monocle3)
  library(dplyr)
  library(ggplot2)
})



base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
# output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA609-v6/monocle/'
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/monocle_v2/'
dir.create(output_dir)
datatag <- 'SA609'
obs_clones <- c('R')
obs_treatment_st <- c('UT','UTT','UTTT','UTTTT','UTTTTT')
obs_clones_untreated <- c('H')
obs_untreated_st <- c('U','UU','UUU','UUUU','UUUUU')

sce <- readRDS(paste0(base_dir,'rnaseq_v6/SA609-v6/SA609_sctransform_normalized.rds'))
# sce <- readRDS(paste0(base_dir,'rnaseq_v6/normalization_evaluation/SA609/SA609_norm_batch_corrected_sce.rds'))
dim(sce)
meta_data <- as.data.frame(colData(sce))
print(table(meta_data$clone, meta_data$treatmentSt))
dim(meta_data)
cells_use_g1 <- meta_data %>%
  dplyr::filter(treatmentSt %in% obs_treatment_st &
                  clone %in% obs_clones) %>%
  dplyr::pull(cell_id)


cells_use_g1 <- meta_data %>%
  dplyr::filter(grepl('T',treatmentSt) &
                  clone %in% obs_clones) %>%
  dplyr::pull(cell_id)
print(length(cells_use_g1))
# cells_use_g2 <- meta_data_tmp %>%
#   dplyr::filter(treatmentSt %in% obs_untreated_st &
#                   clone %in% obs_clones_untreated) %>%
#   dplyr::pull(cell_id)
# print(length(cells_use_g2))

# obs_ts <- c('UTT','UTTT','UTTU')
# obs_ts <- c('UTTT','UTTTT','UTTTU')

cells_use_g1 <- meta_data %>%
  dplyr::filter(treatmentSt %in% obs_ts) %>%
  dplyr::pull(cell_id)

print(length(cells_use_g1))
meta_data <- meta_data %>%
  dplyr::filter(treatmentSt %in% obs_ts)
dim(meta_data)
meta_data <- meta_data %>%
  dplyr::filter(cell_id %in% cells_use_g1)

# grepl('T',treatmentSt)


# rm(tsce)
# tsce <- sce[,c(cells_use_g1,cells_use_g2)]
ext_cells <- downsample_data(meta_data, downsample_ratio=0.8)
# tsce <- sce[,ext_cells]
tsce <- sce[,cells_use_g1]
print(dim(tsce))
table(tsce$treatmentSt, tsce$clone)



# tsce <- readRDS(paste0(base_dir,'rnaseq_v6/SA609-v6/SA609_norm_batch_corrected_sce_Rx.rds'))


if(is.null(rowData(tsce)$Symbol)){
  genes_symb_df <- read.csv(paste0(base_dir,'biodatabase/meta_genes.csv'), check.names = F, stringsAsFactors = F)
  dim(genes_symb_df)
  rownames(genes_symb_df) <- genes_symb_df$gene_ens
  rowData(tsce)$gene_short_name <- genes_symb_df[rownames(tsce),'gene_symb']
}
rowData(tsce)$gene_short_name[1]
tsce$drug_dose <- stringr::str_count(tsce$treatmentSt, "T")
# saveRDS(tsce,file=paste0(base_dir,'rnaseq_v6/SA609-v6/SA609_sctransform_sce_',paste(obs_ts, collapse = '_'),'.rds'))
saveRDS(tsce,file=paste0(base_dir,'rnaseq_v6/SA609-v6/SA609_sctransform_sce_Rx_cloneR.rds'))


samples <- unique(meta_data$sample)
cells_use <- meta_data$Barcode
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA609-v6/'
meta_data <- as.data.frame(colData(tsce))
cds1 <- create_monocle_obj_from_sce(sce_combine, save_dir=NULL, exp='counts', preprocess=F, save_data=F)
summary(cds1@colData$Size_Factor)
meta_cells <- as.data.frame(cds1@colData)
class(meta_cells)
rownames(meta_cells)[1]
# tsce <- tsce[,tsce$treatmentSt!='UT'] 

# rowData(tsce)$gene_id <- rownames(tsce)
# 
# tsce <- tsce[rowData(tsce)$gene_short_name %in% de_genes$gene_symbol,]
# # rowData(sce)$gene_short_name <- rowData(sce)$Symbol
# unique(tsce$clone)
# unique(tsce$treatmentSt)
# print(dim(tsce))
# 


# unique(tsce$drug_dose)
# # rm(sce)
# dim(expression_matrix)
# 
# # de_genes <- readRDS(paste0(output_dir,'SA609_res_monocle_modules_genes.rds')) 
# # de_genes <- readRDS(paste0(output_dir,'SA609_res_monocle.rds')) 
# # class(de_genes)
# # de_genes <- as.matrix(de_genes) %>% as.data.frame()
colnames(cds)[1]
cds <- create_monocle_obj_from_sce(tsce, save_dir=NULL, exp='counts', preprocess=F, save_data=F)
colData(cds)$Size_Factor <- NULL
colData(cds)$Size_Factor <- meta_cells[colnames(cds),'Size_Factor']
saveRDS(cds, paste0(output_dir,'cds_raw_Rx_RxH.rds'))
# cds <- create_monocle_obj_from_sce(tsce, output_dir)
# dim(cds)
# FM <- as.matrix(logcounts(tsce))
# cds <- preprocess_cds(cds, num_dim = 50, norm_method='none')  # data have already processed
cds <- preprocess_cds_v2(cds, as.matrix(logcounts(tsce)), method ='PCA', num_dim=50, scaling = TRUE, verbose=TRUE)
dim(cds)

cds <- reduce_dimension(cds, max_components = 15, preprocess_method = 'PCA',reduction_method = 'UMAP')

saveRDS(cds_backup, paste0(output_dir,datatag,'cds_preprocess_Rx_RxH.rds'))

saveRDS(cds, paste0(output_dir,datatag,'cds_results_Rx.rds'))
# rm(res)
res <- plot_raw_data(cds, output_dir, cls_resolution=0.3)
# cds <- cluster_cells(cds, resolution=0.01, cluster_method="leiden", reduction_method = 'UMAP')
# meta_data <- as.data.frame(colData(sce))
cds <- res$cds
unique(partitions(cds))
cds <- cds[,partitions(cds)!=4]
dim(cds)
meta_cells <- as.data.frame(cds@colData)
unique(meta_cells$treatmentSt)
cds_backup <- cds
cds <- cds[,!grepl('TU',colData(cds)$treatmentSt)]
cds <- cds[,colData(cds)$treatmentSt!='UT']
dim(cds)
# cds <- cds[,colData(cds)[, 'treatmentSt'] == time_bin]
# rowData(cds)
res <- plot_raw_data(cds, output_dir, cls_resolution=0.3)
# t <- colData(cds)
# unique(partitions(res$cds))
cds <- res$cds
dim(cds)
cds <- learn_graph(cds, learn_graph_control=list(ncenter=50,prune_graph=T))
# cds <- cds[,partitions(cds)==1]
# cds <- learn_graph(cds)
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds, time_bin='UTT', feature_use="treatmentSt"))
dim(cds)
p2 <- plot_cells(cds,
                 color_cells_by = "treatmentSt",
                 trajectory_graph_segment_size = 0.75,
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE, 
                 group_label_size = 4, 
                 cell_size = 1)
p2
p22 <- p2
p22 <- p22 + theme(legend.position = "right")
p22
p3 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size=3.5,
                 cell_size = 1)

p3
p <- cowplot::plot_grid(p22, p3, ncol=2, align = 'h', rel_widths = c(1,1))
png(paste0(output_dir, datatag,"_pseudotime_Rx.png"), height =2*350, width= 2*900,res = 2*72)
print(p)
dev.off()
# cds <- readRDS(paste0(output_dir,'cds_Rx.rds'))
# t <- colData(cds)
# unique(t$drug_dose)

de_genes = data.table::fread(paste0(output_dir,'SA609_all_changing_genes2.csv')) %>% as.data.frame()
head(de_genes)
thrsX <- quantile(abs(de_genes$slopeX),0.6)  # X is treated line
thrsY <- quantile(abs(de_genes$slopeY),0.6)  # Y is untreated line
de_genes <- de_genes %>%
  dplyr::filter(abs(slopeX)>=thrsX & abs(slopeY)>=thrsY)

dim(de_genes)
dim(cds)
cds_backup <- cds
dim(cds_backup)
cds <- cds_backup
cds <- cds[,]
cds <- cds[rowData(cds)$gene_short_name %in% de_genes$gene_symbol,]
dim(cds)
gene_fits <- fit_models(cds, model_formula_str = "~drug_dose")

# 
# tsce <- tsce[rownames(tsce) %in% de_x3$ensembl_gene_id,]
# fit_coefs1 <- fit_coefs
# fit_coefs1$gene_short_name
# fit_coefs1 <- fit_coefs1 %>%
#   dplyr::filter(gene_short_name %in% de_genes$gene_symbol)
# dim(fit_coefs1)
res_de <- fit_regression_DE_genes(cds, output_dir, feature_use='drug_dose')
res_de <- res
# create_monocle_obj(expression_matrix, cell_metadata, 
#                    gene_metadata, de_x3,
#                    0.3, output_dir)  




# dim(cds)
# saveRDS(cds, paste0(output_dir,'SA609_monocle_cds_R_Rx_H_UnRx.rds'))
rm(cds)
cds <- readRDS(paste0(output_dir,'SA609_monocle_cds_R_Rx_H_UnRx.rds'))
# de_genes <- res$de_genes
de_genes = data.table::fread(paste0(output_dir,'genes_signif_treated.csv.gz')) %>% as.data.frame()
head(de_genes)
dim(de_genes)
# rowData(cds)$gene_short_name[1]
if(!is.null(de_genes)){
  cds <- cds[rowData(cds)$gene_short_name %in% unique(de_genes$gene_short_name),]  
}else{  # use all genes
  # cds <- cds[1:5,] 
  stop('Load the list of DE genes first')
}
print('DEBUG')
print(dim(cds))
save_dir <- output_dir
find_coregulated_genes(cds, cls_res=0.2, output_dir)

# cds <- cds[,colData(cds)[, 'clone'] =='R']
# dim(cds)
# unique(clusters(cds))
# de_x3 = read.csv(paste0(base_dir,'SA609_rna/deg_analysis/SA609-v6/SA609_UTTT_R_UUUU_H/signif_genes.csv'))
# head(de_x3)

p_pseudo1 <- plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name %in% pr_deg_ids[1:6],],
                                      color_cells_by="pseudotime",
                                      min_expr=0.5)

p_pseudo2 <- plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name %in% pr_deg_ids[1:6],],
                                      color_cells_by="treatmentSt",
                                      min_expr=0.5)

p_total <- cowplot::plot_grid(p_pseudo1, p_pseudo2, ncol=2, align = 'hv')
png(paste0(output_dir, datatag,"_pseudotime_genesplt.png"), height =2*900, width= 2*1000,res = 2*72)
print(p_total)
dev.off()



gene_module_df <- data.table::fread(paste0(output_dir,'SA609_morangene_module_treated.csv.gz')) %>% as.data.frame()
dim(gene_module_df)
gene_module_df$gene_short_name <- genes_symb_df[gene_module_df$id,'gene_symb']
gene_module_df$module[1]
unique(gene_module_df$supermodule)
obs_modules <- c(12, 4, 9, 19, 17, 20, 18, 25, 15, 7, 11, 16, 1, 6)
gene_module_df_backup <- gene_module_df
gene_module_df <- gene_module_df_backup
obs_modules <- c(2,1,4,6)
unique(gene_module_df$module)
obs_modules <- c(36)
gene_module_df$module <- paste0("Module ",gene_module_df$module)
gene_module_df <- gene_module_df %>%
  dplyr::filter(module %in% obs_modules)

gene_module_df <- gene_module_df %>%
  dplyr::filter(module %in% paste0("Module ",obs_modules))
dim(gene_module_df)

p2 <- plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% paste0("Module ",obs_modules)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE,
           cell_size = 0.8)

p1 <- plot_cells(cds,
                 genes=gene_module_df %>% filter(module %in% paste0("Module ",4)),
                 label_cell_groups=FALSE,
                 show_trajectory_graph=FALSE,
                 cell_size = 0.8)

p <- cowplot::plot_grid(p3, p1,p2, ncol=3, rel_widths = c(1,1,1))
png(paste0(output_dir, datatag,"_pseudotime_module",obs_modules,".png"), height =2*250, width= 2*950,res = 2*72)
print(p)
dev.off()

pg <- p4
pg <- plot_cells(cds,
                 genes=pr_deg_ids[1:6],
                 label_cell_groups=FALSE,
                 show_trajectory_graph=FALSE,
                 cell_size = 0.8)
png(paste0(output_dir,"pseudotime_genes.png"), height =2*500, width= 2*850,res = 2*72)
print(pg)
dev.off()

p <- plot_genes_violin(cds[rowData(cds)$gene_short_name %in% pr_deg_ids[1:6],], group_cells_by='treatmentSt', ncol=3, normalize=T, log_scale=T) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
png(paste0(output_dir, "genes_violin_pseudo.png"), height =2*500, width= 2*750,res = 2*72)
print(p)
dev.off()


gene_module_df <- gene_module_df %>% inner_join(genome_genes_df, by=c('id'='Ensembl'))
obs_genes <- unique(gene_module_df$Symbol)
length(obs_genes)
obs_genes[1:6]
save_dir <- output_dir

paste(obs_genes, collapse = ', ')
ref_genes <- c()
for(pw in names(pathway_sets)){
  ref_genes <- c(ref_genes, unlist(pathway_sets[[pw]]))
}  
ref_genes <- unique(ref_genes)
sum(obs_genes %in% ref_genes)


cds_sub <- choose_graph_segments(cds)
# plot_cells(cds_sub,
#            color_cells_by = "pseudotime",
#            label_cell_groups=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE,
#            graph_label_size=1.5)
cds_sub <- order_cells(cds_sub, root_pr_nodes=get_earliest_principal_node(cds, time_bin='UTT', feature_use="treatmentSt"))

cds_sub <- cds_sub[gene_module_df$id,]
dim(cds_sub)
subset_pr_test_res <- graph_test(cds_sub, neighbor_graph="principal_graph", cores=4)

plot_genes_in_pseudotime(cds_sub,
                         color_cells_by="drug_dose",
                         min_expr=0.5)

rm(cds_subset)
cds_subset <- choose_cells(cds)
cds_subset <- choose_graph_segments(cds)
dim(cds_subset)
subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
subset_pr_test_res <- subset(subset_pr_test_res, q_value < 0.05)
subset_pr_test_res <- subset_pr_test_res[order(subset_pr_test_res$sct.variance, decreasing=T),]
subset_pr_test_res$gene_short_name[1:4]

subset_pr_test_res <- subset_pr_test_res %>%
  dplyr::filter(subset_pr_test_res$gene_short_name %in% gene_module_df$gene_short_name)

pr_deg_ids <- row.names(subset_pr_test_res)
length(subset_pr_test_res$gene_short_name)

gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.001)

agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, 
                                levels = row.names(agg_mat)[module_dendro$order])

plot_cells(cds_subset,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

p1 <- plot_cells(cds_subset, genes=subset_pr_test_res$gene_short_name[2:4],
           show_trajectory_graph=F,
           label_cell_groups=T,
           label_leaves=FALSE)

p2 <- plot_cells(cds_subset,
                 color_cells_by = "treatmentSt",
                 label_groups_by_cluster=T,
                 show_trajectory_graph=T,
                 label_leaves=FALSE,
                 label_branch_points=FALSE)

p <- cowplot::plot_grid(p2,p1, ncol=2, rel_widths = c(1,3))
png(paste0(output_dir,datatag,"_pseudotime_",paste(obs_ts, collapse = '_'),".png"), height =2*350, width= 2*800,res = 2*72)
print(p)
dev.off()
t1 <- cds_subset[rowData(cds_subset)$gene_short_name %in% subset_pr_test_res$gene_short_name[2:4],]
dim(t1)
colData(t1)$Size_Factor <- 1
plot_genes_in_pseudotime(t1,
                         color_cells_by="treatmentSt")
