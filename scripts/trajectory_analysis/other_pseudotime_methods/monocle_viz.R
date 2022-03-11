suppressPackageStartupMessages({
  library(monocle3)
  library(dplyr)
  library(ggplot2)
})
get_earliest_principal_node <- function(cds, time_bin="UT", feature_use="treatmentSt"){
  cell_ids <- which(colData(cds)[, feature_use] == time_bin)
  
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
create_monocle_obj_v2 <- function(expression_matrix, cell_metadata, 
                                  gene_metadata, save_dir=NULL){
  cds <- new_cell_data_set(expression_matrix,
                           cell_metadata=cell_metadata,
                           gene_metadata = gene_metadata)
  print(dim(cds))
  
  cds <- preprocess_cds(cds, num_dim = 50, norm_method='none')  # data have already processed
  cds <- reduce_dimension(cds, preprocess_method = 'PCA')
  # print('Save data into folder...')
  # saveRDS(cds, paste0(save_dir,'cds_Rx.rds'))
  return(cds)
}  
create_monocle_obj_from_sce <- function(tsce, save_dir=NULL){
  expression_matrix <- as.matrix(logcounts(tsce))
  cell_metadata <- as.data.frame(colData(tsce))
  gene_metadata <- as.data.frame(rowData(tsce))
  cds <- new_cell_data_set(expression_matrix,
                           cell_metadata=cell_metadata,
                           gene_metadata = gene_metadata)
  print(dim(cds))
  
  cds <- preprocess_cds(cds, num_dim = 50, norm_method='none')  # data have already processed
  cds <- reduce_dimension(cds, preprocess_method = 'PCA')
  # print('Save data into folder...')
  # saveRDS(cds, paste0(save_dir,'cds_Rx.rds'))
  return(cds)
}  
plot_raw_data <- function(cds, save_dir, cls_resolution=0.3){
  cds <- cluster_cells(cds, resolution=cls_resolution, cluster_method="leiden")
  p1 <- plot_cells(cds, color_cells_by = "partition",
                   label_groups_by_cluster=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE)
  
  # cds <- learn_graph(cds)
  p2 <- plot_cells(cds,
                   color_cells_by = "treatmentSt",
                   label_groups_by_cluster=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE)
  
  p <- cowplot::plot_grid(p1, p2, ncol=2, align = 'h') #label=c('Drug Cycle Treatment','Untreated')
  png(paste0(save_dir, "umap_inputdata.png"), height =2*400, width= 2*900,res = 2*72)
  print(p)
  dev.off()
  
  return(list(p=p, cds=cds))
}
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA609-v6/monocle/'
# cds <- readRDS(paste0(output_dir,'cds_Rx.rds'))
# res <- plot_raw_data(cds, output_dir, cls_resolution=0.3)

datatag <- 'SA609'
obs_clones <- c('R')
# obs_treatment_st <- c('UT','UTT','UTTT','UTTTT','UTTTTT')

sce <- readRDS(paste0(base_dir,'rnaseq_v6/normalization_evaluation/SA609/SA609_norm_batch_corrected_sce.rds'))
dim(sce)
meta_data_tmp <- as.data.frame(colData(sce))
print(table(meta_data_tmp$clone, meta_data_tmp$treatmentSt))
# cells_use_g1 <- meta_data_tmp %>%
#   dplyr::filter(treatmentSt %in% obs_treatment_st &
#                   clone %in% obs_clones) %>%
#   dplyr::pull(cell_id)

cells_use_g1 <- meta_data_tmp %>%
  dplyr::filter(clone %in% obs_clones) %>%
  dplyr::pull(cell_id)

tsce <- sce[,cells_use_g1]
print(dim(tsce))

if(is.null(rowData(tsce)$Symbol)){
  genes_symb_df <- read.csv(paste0(base_dir,'biodatabase/meta_genes.csv'), check.names = F, stringsAsFactors = F)
  dim(genes_symb_df)
  rownames(genes_symb_df) <- genes_symb_df$gene_ens
  rowData(tsce)$gene_short_name <- genes_symb_df[rownames(tsce),'gene_symb']
}
rowData(tsce)$gene_short_name[1]
tsce$drug_dose <- 2 * stringr::str_count(tsce$treatmentSt, "T")
# unique(tsce$drug_dose)

expression_matrix <- as.matrix(logcounts(tsce))
cell_metadata <- as.data.frame(colData(tsce))
gene_metadata <- as.data.frame(rowData(tsce))

cds <- create_monocle_obj_v2(expression_matrix, cell_metadata, gene_metadata, output_dir)
cds <- cluster_cells(cds, resolution=0.3, cluster_method="leiden")
cds <- learn_graph(cds)
saveRDS(cds, paste0(output_dir,'cds_cloneR.rds'))
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds, time_bin="UT", feature_use="treatmentSt"))

p1 <- plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

p2 <- plot_cells(cds,
                 color_cells_by = "treatmentSt",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE)

p <- cowplot::plot_grid(p1, p2, ncol=2, align = 'h') #label=c('Drug Cycle Treatment','Untreated')
png(paste0(save_dir, "umap_trajectory_drug_holiday.png"), height =2*400, width= 2*900,res = 2*72)
print(p)
dev.off()

saveRDS(cds, paste0(output_dir,'cds_cloneR.rds'))