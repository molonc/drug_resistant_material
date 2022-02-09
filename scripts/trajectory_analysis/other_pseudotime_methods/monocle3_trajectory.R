
library(monocle3)
library(dplyr)
library(ggplot2)

expression_matrix <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_expression.rds"))
cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_colData.rds"))
gene_annotation <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_rowData.rds"))

# dim(expression_matrix)
# dim(cell_metadata)
# ncells <- 100
# exp <- log2(expression_matrix[,1:ncells] + 1)
# cds <- new_cell_data_set(as.matrix(exp),
#                          cell_metadata = cell_metadata[1:ncells,],
#                          gene_metadata = gene_annotation)
# 

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
dim(cds)
t <- colData(cds)
unique(t$embryo.time)
fit_coefs <- coefficient_table(gene_fits)
View(head(fit_coefs))

# expression_matrix[1:3,1:3]
class(expression_matrix)
class(cell_metadata)
class(gene_annotation)
# cds <- preprocess_cds(cds, num_dim = 20, norm_method='none')
cds <- preprocess_cds(cds, num_dim = 50)
# cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")

cds <- reduce_dimension(cds, preprocess_method = 'PCA')
# plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cell.type")
unique(t$embryo.time.bin)
ciliated_genes <- c("che-1",
                    "hlh-17",
                    "nhr-6",
                    "dmd-6",
                    "ceh-36",
                    "ham-1")
cds_subset <- cds[rowData(cds)$gene_short_name %in% ciliated_genes,]
gene_fits <- fit_models(cds_subset, model_formula_str = "~embryo.time")
gene_fits2 <- fit_models(cds_subset, model_formula_str = "~embryo.time.bin")
# View(gene_fits2[1:10,])
# gene_fits$term
# plot_cells(cds,
#            genes=ciliated_genes,
#            label_cell_groups=FALSE,
#            show_trajectory_graph=FALSE)

cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

cds <- learn_graph(cds)
# plot_cells(cds,
#            color_cells_by = "cell.type",
#            label_groups_by_cluster=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE)

# plot_cells(cds,
#            color_cells_by = "embryo.time.bin",
#            label_cell_groups=FALSE,
#            label_leaves=TRUE,
#            label_branch_points=TRUE,
#            graph_label_size=1.5)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="130-170"){
  cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

cds_subset <- cds[rowData(cds)$gene_short_name %in% ciliated_genes,]


# gene_fits <- fit_models(cds_subset, model_formula_str = "~embryo.time + batch")
# fit_coefs <- coefficient_table(gene_fits)
# fit_coefs %>% filter(term != "(Intercept)") %>%
#   select(gene_short_name, term, q_value, estimate)

rownames(cds_subset)[1]
t <- colData(cds)
class(t)
head(t[,1:4])
colnames(t)
cds_subset <- cds[,metadata(cds)]


# Genes modules
pr_graph_test_res <- graph_test(cds, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
length(pr_deg_ids)
# gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(0,10^seq(-6,-1)))
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-2)


cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=partitions(cds)) #[colnames(cds)]
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)

# Genes change as function of pseudotime with neighbor_graph="principal_graph"
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))


gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(0,0.01))

plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% c(27, 10, 7, 30)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

AFD_genes <- c("gcy-8", "dac-1", "oig-8")
AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,
                       colData(cds)$cell.type %in% c("AFD")]

plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="embryo.time.bin",
                         min_expr=0.5)


# evaluate_fits(gene_fits)
# time_batch_models <- fit_models(cds_subset,
#                                 model_formula_str = "~embryo.time + batch",
#                                 expression_family="negbinomial")
# time_models <- fit_models(cds_subset,
#                           model_formula_str = "~embryo.time",
#                           expression_family="negbinomial")
# compare_models(time_batch_models, time_models) %>% select(gene_short_name, q_value)




# pt <- partitions(cds)
# class(pt)
# unique(pt)

# cds_p1 <- cds[,partitions(cds)=='1']
# dim(cds_p1)
# colData(cds_p1)
t <- summary(as.factor(colData(cds_p1)[, feature_use]))
cond <- colData(cds_p1)[, feature_use] %in% c('UTTTU','UTTTT')
cds_p1 <- cds_p1[,!cond]
cds_p1 <- cluster_cells(cds_p1)
plot_cells(cds_p1, color_cells_by = "partition",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
cds_p1 <- learn_graph(cds_p1)
plot_cells(cds_p1,
           color_cells_by = "treatmentSt",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)


cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
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
# p <- p1 + p2
# cds_sub1 <- choose_graph_segments(cds_p1)
# cds_sub1 <- choose_cells(cds_p1)
# dim(cds_sub1)
summary(as.factor(partitions(cds_sub1)))
# pr_graph_test_res <- graph_test(cds_sub1, neighbor_graph="knn", cores=8)
# cds_sub1 <- learn_graph(cds_sub1)
subset_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05 & p_value <0.05))
# cds_pr_test_res <- graph_test(cds_sub1, neighbor_graph="principal_graph", cores=4)
# pr_deg_ids <- row.names(subset(cds_pr_test_res, q_value < 0.05))
length(pr_deg_ids)
subset_pr_test_res$p_value[1]
saveRDS(agg_mat, paste0(base_dir,'rnaseq_v6/SA609-v6/SA609_res_monocle_modules_genes.rds'))
unique(colnames(agg_mat))
cds_p2 <- cluster_cells(cds_p1)
summary(as.factor(partitions(cds_p2)))
?cluster_cells()
clusters(cds_p1)[1:3]
View(head(agg_mat))

plot_cells(cds_sub1,
           genes=gene_module_df %>% filter(module %in% c(53,46,52)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


t <- colData(cds_subset) %>% as.data.frame()
unique(t$treatmentSt)
cds_subset <- cds_p1