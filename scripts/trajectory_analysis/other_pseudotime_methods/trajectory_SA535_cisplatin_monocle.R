


suppressPackageStartupMessages({
  library(monocle3)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
})

data.table::fread()
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA535-v6/'
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/monocle_v2/'
dir.create(output_dir)
datatag <- 'SA535'
# obs_clones <- c('R')
obs_treatment_st <- c('UT','UTT','UTTT','UTTTT','UTTTTT')
# obs_clones_untreated <- c('H')
obs_untreated_st <- c('U','UU','UUU','UUUU','UUUUU')

meta_cells <- data.table::fread('/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/snakemake/metasample_SA535_cisplatin_untreated.csv') %>% as.data.frame()
View(meta_cells)
meta_cells$treatmentSt
meta_cells <- meta_cells %>%
  dplyr::filter(PDX!='SA535_untreated')
meta_cells <- meta_cells %>%
  dplyr::select(mouse_id, batch_info, passage)
meta_cells <- meta_cells %>%
  dplyr::filter(!grepl('TU',treatmentSt))
rownames(meta_cells) <- meta_cells$library_id
sce_combine <- load_data(meta_cells, base_dir,tag='CIS_Rx_RxH')
sce_combine$treatmentSt <- NULL
sce_combine$treatmentSt <- meta_cells[sce_combine$library_id,'treatmentSt']
unique(sce_combine$library_id) %in% meta_cells$library_id

sce <- sce_combine[,sce_combine$treatmentSt %in% unique(meta_cells$treatmentSt)]
unique(sce$treatmentSt)
# sce <- readRDS(paste0(base_dir,'SA535_CIS_CX_untreated_total_sce_v3.rds'))
# sce <- readRDS(paste0(base_dir,'rnaseq_v6/normalization_evaluation/SA609/SA609_norm_batch_corrected_sce.rds'))
saveRDS(sce, file=paste0(base_dir,"raw_SA535_CIS_Rx.rds"))
sce_raw <- readRDS(paste0(base_dir,"scTransform_SA535_CIS_Rx.rds"))
dim(sce_raw)
sum(colnames(cds) %in% colnames(sce_raw))
rownames(cds)[1]
rownames(sce_raw)[1]
# rownames(sce_raw) <- rowData(sce_raw)$ID

norm_mtx <- logcounts(sce_raw[rownames(cds),colnames(cds)])


output_file <- paste0(base_dir,"scTransform_SA535_CIS_Rx.rds")
zero_cbind <- DelayedArray::rowMeans(assay(sce_combine, 'counts') == 0)
cut_off_overall <- 0.05
sce_combine <- sce_combine[names(zero_cbind[zero_cbind <= (1 - cut_off_overall)]), ]
dim(sce_combine)
rownames(sce_combine) <- rowData(sce_combine)$ID
rownames(sce_combine)[1]
colnames(sce_combine)[1]

sc_norm <- normalize_SCTransform(sce_combine, base_dir, paste0(datatag,'_Rx_RxH'), return_data=T, output_file)
sc_norm <- readRDS(paste0(base_dir,'scTransform_SA535_CIS_Rx.rds'))
norm_sce <- sc_norm[,colnames(sce)]
dim(sce)
dim(norm_sce)
rownames(sce) <- rowData(sce)$ID
colnames(sce)[1]
tsce <- sce
rowData(tsce)$gene_short_name <- rowData(tsce)$Symbol
expression_matrix <- as.matrix(counts(tsce))
cell_metadata <- as.data.frame(colData(tsce))
gene_metadata <- as.data.frame(rowData(tsce))
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata=cell_metadata,
                         gene_metadata = gene_metadata)
print(dim(cds))
class(counts(tsce))
FM <- as.matrix(logcounts(norm_sce))
# cds <- preprocess_cds(cds, num_dim = 50, norm_method='none')  # data have already processed
cds <- preprocess_cds_v2(cds, as.matrix(logcounts(norm_sce)), method ='PCA', num_dim=50, scaling = TRUE, verbose=TRUE)
dim(cds)

cds <- reduce_dimension(cds, max_components = 15, preprocess_method = 'PCA',reduction_method = 'UMAP')

saveRDS(cds, paste0(output_dir,datatag,'cds_Rx_raw.rds'))
res <- plot_raw_data(cds, output_dir, cls_resolution=0.3)
cds <- cluster_cells(cds, resolution=0.005, cluster_method="leiden", reduction_method = 'UMAP')
meta_data <- as.data.frame(colData(sce))
cds <- cluster_cells(cds, resolution=0.001, cluster_method="leiden", reduction_method = 'UMAP')
unique(clusters(cds))
unique(colData(cds)$cluster_label)
colData(cds)$cluster_label <- paste0('Cls_',clusters(cds))
print(table(meta_data$clone, meta_data$treatmentSt))
dim(meta_data)
cells_use_g1 <- meta_data %>%
  dplyr::filter(treatmentSt %in% obs_treatment_st &
                  clone %in% obs_clones) %>%
  dplyr::pull(cell_id)
print(length(cells_use_g1))
plot_cells(cds, color_cells_by = "partition")
cds <- learn_graph(cds)
# cds <- learn_graph(cds)
p1 <- plot_cells(cds,
                 # reduction_method='PCA',
           color_cells_by = "treatmentSt",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

p1 <- plot_cells(cds,
                 color_cells_by = "cluster_label",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 group_label_size = 4, 
                 cell_size = 1)

p1

extract_DE_genes_clusters <- function(cds, genes_df, output_dir, cls_ls, feature_use='cluster_label')
cds_backup <- cds

cds <- cds[,colData(cds)$treatmentSt!='UT']
cds <- learn_graph(cds, learn_graph_control=list(ncenter=200))
# saveRDS(cds, paste0(output_dir,datatag,'cds_Rx_processed.rds'))
saveRDS(cds, paste0(output_dir,datatag,'cds_Rx_processed_v2.rds'))
cds <- readRDS(paste0(output_dir,datatag,'_cds_Rx_processed_v2.rds'))
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
png(paste0(output_dir, datatag,"_pseudotime.png"), height =2*350, width= 2*900,res = 2*72)
print(p)
dev.off()


p <- cowplot::plot_grid(p1, p22, p3, ncol=3, align = 'h', rel_widths = c(1,1.2,1.2))
png(paste0(output_dir, datatag,"_pseudotime_clusters.png"), height =2*300, width= 2*1000,res = 2*72)
print(p)
dev.off()
dim(genes_df)
colData(cds)$cluster_label <- as.factor(colData(cds)$cluster_label)
cls_ls <- c('Cls_1','Cls_3')
de13 = data.table::fread(paste0(output_dir,paste(cls_ls, collapse='_'),'_DE_genes.csv.gz')) %>% as.data.frame()
extract_DE_genes_clusters(cds, genes_df, output_dir, cls_ls, feature_use='cluster_label')
cls_ls <- c('Cls_1','Cls_2')
de12 = data.table::fread(paste0(output_dir,paste(cls_ls, collapse='_'),'_DE_genes.csv.gz')) %>% as.data.frame()
extract_DE_genes_clusters(cds, genes_df, output_dir, cls_ls, feature_use='cluster_label')

cls_ls <- c('Cls_2','Cls_5')
de25 = data.table::fread(paste0(output_dir,paste(cls_ls, collapse='_'),'_DE_genes.csv.gz')) %>% as.data.frame()
extract_DE_genes_clusters(cds, genes_df, output_dir, cls_ls, feature_use='cluster_label')
colData(cds)$cluster_label <- ifelse(colData(cds)$cluster_label %in% c('Cls_2','Cls_5'),'Cls_25',colData(cds)$cluster_label)
cls_ls <- c('Cls_1','Cls_25')
de125 = data.table::fread(paste0(output_dir,paste(cls_ls, collapse='_'),'_DE_genes.csv.gz')) %>% as.data.frame()
extract_DE_genes_clusters(cds, genes_df, output_dir, cls_ls, feature_use='cluster_label')
dim(cds)

genes_cls <- union(de13$gene_id, de12$gene_id)
genes_cls <- union(genes_cls, de25$gene_id)
length(genes_cls)
genes_cls <- intersect(unique(genes_cls), genes_df$gene_id)
genes_use <- genes_cls
datatag <- paste0(datatag,'_DE_genes_clusters')
find_coregulated_genes_v2(cds, genes_cls, output_dir, paste0(datatag,'_DE_genes_clusters'), cls_res=0.2)
de12 <- de12[order(abs(de12$estimate), decreasing=T),]

unique(colData(cds)$cluster_label)

  
cds_subset <- choose_cells(cds)
cds_subset <- cds_subset[rownames(cds_subset) %in% de_genes$gene_id,]
dim(cds_subset)

# cds_subset <- choose_graph_segments(cds)
dim(cds_subset)
# unique(colData(cds_subset)$treatmentSt)
cds_pr_res <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
# cds_pr_res <- graph_test(cds_subset, cores=5)
cds_pr_res <- subset(cds_pr_res, q_value < 0.05)
dim(cds_pr_res)
dim(cds)
View(head(cds_pr_res))
summary(cds_pr_res$morans_I)
dim(cds)
cds_pr_res <- cds_pr_res[order(cds_pr_res$morans_I, decreasing=T),]

cds_pr_res1 <- cds_pr_res[cds_pr_res$morans_I>=quantile(cds_pr_res$morans_I,0.75),]
dim(cds_pr_res1)
data.table::fwrite(cds_pr_res, paste0(output_dir,'pseudotime_genes_.csv.gz'))
cds_pr_res <- data.table::fread(paste0(output_dir,'pseudotime_genes_.csv.gz')) %>% as.data.frame()
View(cds_pr_res[1:10,])
colnames(cds_pr_res)
dim(cds_pr_res)
cds_pr_res <- cds_pr_res %>%
     dplyr::select(p_value,Symbol, morans_I,biotype, description)
data.table::fwrite(cds_pr_res1, paste0(output_dir,'pseudotime_genes_thrs_0.75qt.csv.gz'))
cds_pr_res1 <- cds_pr_res1[order(cds_pr_res1$morans_I, decreasing=T),]


pr_deg_ids <- cds_pr_res1$gene_short_name
pr_deg_ids <- pr_deg_ids[pr_deg_ids!='MALAT1']
p4 <- plot_cells(cds, genes=pr_deg_ids[1:6],
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
p4

genes_df <- get_expr_ratio_across_time(cds, exprs='counts', min_pct=0.1)
dim(genes_df)
head(genes_df)
genes_df$exp_quality
genes_df <- genes_df %>%
  dplyr::filter(exp_quality==T)
summary(genes_df$pct_expr)

find_coregulated_genes_v2(cds, de_genes$gene_id, output_dir, paste0(datatag,'_moran'), cls_res=0.2)

module_genes = data.table::fread(paste0(output_dir,datatag,'gene_module_treated.csv.gz')) %>% as.data.frame()
dim(module_genes)
module_genes = data.table::fread(paste0(output_dir,'SA535_v2gene_module_treated.csv.gz')) %>% as.data.frame()
head(module_genes)
obs_module <- 4
m_genes <- module_genes %>%
  dplyr::filter(module==obs_module)%>%
  dplyr::pull(id)
length(m_genes)
dim(de_genes)
obs_genes <- de_genes %>%
  dplyr::filter(gene_id %in% m_genes)
obs_genes <- obs_genes[order(obs_genes$estimate, decreasing=T),]
dim(obs_genes)
View(de12[1:10,])
de12 <- de12[order(de12$estimate, decreasing=T),]
up_genes <- de12$gene_short_name[1:3]
up_genes <- de12$gene_short_name
de12 <- de12[order(-de12$estimate, decreasing=T),]
down_genes <- de12$gene_short_name[1:3]

p_pseudo1 <- plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name %in% c(up_genes, down_genes),],
                               color_cells_by="pseudotime",
                               min_expr=0.5)

p_pseudo2 <- plot_genes_in_pseudotime(cds[rowData(cds)$gene_short_name %in% c(up_genes, down_genes),],
                                     color_cells_by="treatmentSt",
                                     min_expr=0.5)

p_total <- cowplot::plot_grid(p_pseudo1, p_pseudo2, ncol=2, align = 'hv')
png(paste0(output_dir, datatag,"_pseudotime_cls1_vs_cls3_pseudoplot.png"), height =2*900, width= 2*1000,res = 2*72)
print(p_total)
dev.off()
p <- plot_cells(cds,
           genes=up_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE,
           cell_size = 0.8)
# p
pp <- plot_cells(cds,
                color_cells_by = "pseudotime",
                label_cell_groups=F,
                show_trajectory_graph=T,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                graph_label_size=3.5,
                cell_size = 0.8)

p1 <- plot_cells(cds,
                genes=down_genes,
                label_cell_groups=FALSE,
                show_trajectory_graph=FALSE,
                cell_size = 0.8)
# p1
p_total <- cowplot::plot_grid(p, p1, ncol=1, align = 'v', labels = c('Up','Down'))
png(paste0(output_dir, datatag,"_pseudotime_cls1_vs_cls3.png"), height =2*450, width= 2*800,res = 2*72)
print(p_total)
dev.off()

unique(colData(cds)$clone)


p1 <- plot_genes_in_pseudotime(cds[as.character(obs_genes$gene_id[1:6]),],
                         color_cells_by="pseudotime",
                         min_expr=0.5)
png(paste0(output_dir, datatag,"_pseudotime_module_",obs_module,"_clusters_pseudotime.png"), height =2*900, width= 2*500,res = 2*72)
print(p1)
dev.off()

plot_cells(cds_subset,
           genes=gene_module_df %>% filter(module %in% c(3, 6)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

dim(gene_module_df)

colData(cds)$drug_dose <- stringr::str_count(colData(cds)$treatmentSt, "T")
# filtered_cds <- cds[]
dim(cds)
res <- fit_regression_DE_genes(cds, genes_df, output_dir, feature_use='drug_dose')
# TO DO: num_cells_expressed
de_genes <- res$de_genes
# https://cole-trapnell-lab.github.io/monocle3/docs/differential/
de_genes = data.table::fread(paste0(output_dir,'genes_signif_treated.csv.gz')) %>% as.data.frame()
head(de_genes)
de_genes$gene_id
dim(de_genes)
de_genes$gene_id[1]
de_genes <- de_genes %>%
  dplyr::filter(gene_id %in% genes_df$gene_id)

# cds_subset <- choose_cells(cds)
length(res$de_genes$gene_id)
rownames(cds)[1]
dim(cds)
res$de_genes$term[1]
cds_subset <- cds[res$de_genes$gene_id,]
subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.001)
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, 
                                levels = row.names(agg_mat)[module_dendro$order])

plot_cells(cds_subset,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)