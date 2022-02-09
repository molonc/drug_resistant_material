suppressPackageStartupMessages({
  library(monocle3)
  library(dplyr)
  library(ggplot2)
})
source('/home/htran/Projects/farhia_project/rnaseq/method_testing/monocle_utils.R')

base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA609-v6/monocle/'
datatag <- 'SA609'

obs_ts <- c('UTTT','UTTTT','UTTTU')
obs_ts <- c('UTT','UTTT','UTTTT')
tsce <- readRDS(paste0(base_dir,'rnaseq_v6/SA609-v6/SA609_norm_batch_corrected_sce_',paste(obs_ts, collapse = '_'),'.rds'))
cds <- create_monocle_obj_from_sce(tsce, output_dir)
dim(cds)
# rowData(cds)
res <- plot_raw_data(cds, output_dir, cls_resolution=0.3)
# t <- colData(cds)
# unique(partitions(res$cds))
cds <- res$cds
# cds <- cds[,partitions(cds)==1]
cds <- learn_graph(cds)
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



