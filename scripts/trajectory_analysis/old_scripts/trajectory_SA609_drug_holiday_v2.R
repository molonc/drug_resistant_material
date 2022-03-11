suppressPackageStartupMessages({
  library(monocle3)
  library(dplyr)
  library(ggplot2)
})
source('/home/htran/Projects/farhia_project/rnaseq/method_testing/monocle_utils.R')

base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA609-v6/monocle/'
datatag <- 'SA609'

# obs_ts <- c('UTTT','UTTTT','UTTTU')
obs_ts <- c('UTT','UTTT','UTTU')
# tsce <- readRDS(paste0(base_dir,'rnaseq_v6/SA609-v6/SA609_norm_batch_corrected_sce_',paste(obs_ts, collapse = '_'),'.rds'))

# tsce <- readRDS(paste0(base_dir,'rnaseq_v6/SA609-v6/SA609_sctransform_sce_',paste(obs_ts, collapse = '_'),'.rds'))
# tsce <- readRDS(paste0(base_dir,'rnaseq_v6/SA609-v6/SA609_sctransform_sce_UTT_UTTT_UTTU.rds'))

# meta_data <- as.data.frame(colData(tsce))
# meta_data <- meta_data %>%
#   dplyr::filter(treatmentSt %in% obs_ts)
# dim(meta_data)
# meta_data <- meta_data %>%
#   dplyr::filter(cell_id %in% cells_use_g1)
# ext_cells <- downsample_data(meta_data, downsample_ratio=0.8)

# cds <- create_monocle_obj_from_sce(tsce, output_dir)

cds <- readRDS(paste0(base_dir,'rnaseq_v6/SA609-v6/monocle/SA609_monocle_UTT_UTTT_UTTTT_cds.rds'))
dim(cds)

# rowData(cds)
res <- plot_raw_data(cds, output_dir, cls_resolution=0.3)

cds <- cluster_cells(cds, resolution=0.001, cluster_method="leiden")
unique(clusters(cds))
p1 <- plot_cells(cds, color_cells_by = "partition",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE)

# t <- colData(cds)
# unique(partitions(res$cds))
cds <- res$cds
# cds <- cds[,partitions(cds)==1]
cds <- learn_graph(cds, learn_graph_control=list(ncenter=50))
p1 <- plot_cells(cds,
                 color_cells_by = "treatmentSt",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE)
obs_ts
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds, time_bin=obs_ts[1], feature_use="treatmentSt"))

p2 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size=1.5)
p <- cowplot::plot_grid(p1, p2, ncol=2, align = 'h')
png(paste0(output_dir,datatag,"_pseudotime_",paste(obs_ts, collapse = '_'),".png"), height =2*350, width= 2*800,res = 2*72)
print(p)
dev.off()

saveRDS(cds, paste0(output_dir,datatag,'_monocle_',paste(obs_ts, collapse = '_'),'_cds.rds'))



