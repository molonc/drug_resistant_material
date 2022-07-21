

source('/home/htran/Projects/farhia_project/rnaseq/drug_manuscript/viz_umap_figs/viz_umaps.R')
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/'
output_dir <- paste0(base_dir,'umap_figs/figs/')

# Clone color definition
color_df <- data.table::fread(paste0(base_dir,'umap_figs/colorcode_total_v3.csv'))
# dim(color_df)
# summary(as.factor(color_df$datatag))

# clone_id, timepoint, treatmentstr
## Testing


# c('UnRx','Rx','RxH')
tag <- 'SA609'
library_grouping_fn <- paste0(base_dir,'cell_clones/',tag, '_DLP_library_groupings.csv')
cell_clones_fn <- paste0(base_dir,'cell_clones/',tag, '_cell_clones_total.csv')
cell_clones <- get_meta_data(cell_clones_fn, library_grouping_fn, tag)
# data.table::fwrite(cell_clones, paste0(base_dir,'cell_clones/',tag, '_cell_clones_metadata.csv'))

# Look at clone E
df <- color_df %>%
  dplyr::filter(datatag==tag)
cols_use <- df$colour
names(cols_use) <- df$clone_id
unique(cell_clones$treatment_desc)
res_SA609 <- plot_fill_barplot_wholedataset(cell_clones, cols_use, output_dir, 
                                            datatag, plottitle=NULL, plotlegend=F)
res_SA609$p
res_SA609$plg


tag <- 'SA535'
library_grouping_fn <- paste0(base_dir,'cell_clones/',tag, '_DLP_library_groupings.csv')
cell_clones_fn <- paste0(base_dir,'cell_clones/',tag, '_cell_clones.csv')
cell_clones <- get_meta_data(cell_clones_fn, library_grouping_fn, tag)
# data.table::fwrite(cell_clones, paste0(base_dir,'cell_clones/',tag, '_cell_clones_metadata.csv'))
df <- color_df %>%
  dplyr::filter(datatag==tag)
cols_use <- df$colour
names(cols_use) <- df$clone_id
res_SA535 <- plot_fill_barplot_wholedataset(cell_clones, cols_use, output_dir, 
                                            datatag, plottitle=NULL, plotlegend=F)
# res_SA535$p
# res_SA535$plg

tag <- 'SA1035'
library_grouping_fn <- paste0(base_dir,'cell_clones/',tag, '_DLP_library_groupings.csv')
cell_clones_fn <- paste0(base_dir,'cell_clones/',tag, '_cell_clones.csv')
cell_clones <- get_meta_data(cell_clones_fn, library_grouping_fn, tag)
# data.table::fwrite(cell_clones, paste0(base_dir,'cell_clones/',tag, '_cell_clones_metadata.csv'))
df <- color_df %>%
  dplyr::filter(datatag==tag)
cols_use <- df$colour
names(cols_use) <- df$clone_id
res_SA1035 <- plot_fill_barplot_wholedataset(cell_clones, cols_use, output_dir, 
                                             datatag, plottitle=NULL, plotlegend=F)
# res_SA1035$p
# res_SA1035$plg



## Loading DLP tree
input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA501_v2/tree_viz_dream/'
dlp_501 <- readRDS(paste0(input_dir,'summary_tree.rds'))

input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA530_v2/tree_viz_dream/'
dlp_530 <- readRDS(paste0(input_dir,'summary_tree.rds'))

input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA604_v2/tree_viz_dream/'
dlp_604 <- readRDS(paste0(input_dir,'summary_tree.rds'))

input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_Tyler_v2/tree_viz_dream/'
dlp_1035 <- readRDS(paste0(input_dir,'summary_tree.rds'))

input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_fitness/tree_viz_dream/'
dlp_535 <- readRDS(paste0(input_dir,'summary_tree.rds'))

input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA609/tree_viz_dream/'
dlp_609 <- readRDS(paste0(input_dir,'summary_tree.rds'))


pUnRx <- cowplot::plot_grid(dlp_501,dlp_530, dlp_604, ncol=3, rel_widths = c(1,0.8,1),
                            labels = c('Pt1','Pt2','Pt3'), label_y = 0.2,label_x=0.8, hjust = -0.5)
pSeries <- cowplot::plot_grid(dlp_609,dlp_535, dlp_1035, ncol=3, rel_widths = c(0.8,1,1),
                              labels = c('Pt4','Pt5','Pt6'), label_y = 0.2,label_x=0.75, hjust = -0.5)

pRxDLP <- cowplot::plot_grid(res_SA609$p, res_SA535$p, res_SA1035$p, ncol=3)#rel_heights = c(1,1,1,0.1,0.1,0.1)

# ptotal_DLP <- cowplot::plot_grid(NULL, pUnRx,NULL, pSeries, pRxDLP, nrow=5, rel_heights = c(0.1,1,0.1,1,2),
#                                  labels = c('     -Rx PDX tumors','',
#                                             '     -Rx/+Rx time series PDX tumors',''),
#                                  # , label_x = 0, label_y = 1
#                                  label_x = .01, hjust = 0, label_size=13)+
#   theme(plot.background = element_rect(fill = "white", colour = "white"))

pDLP_tree1 <- cowplot::plot_grid(NULL, pUnRx,NULL, pSeries, nrow=4, rel_heights = c(0.1,1,0.1,1),
                                 labels = c('     UnRx PDX tumors','',
                                            '     UnRx/Rx,RxH time series PDX tumors',''),
                                 # , label_x = 0, label_y = 1
                                 label_x = .01, hjust = 0, label_size=12)+
  theme(plot.background = element_rect(fill = "white", colour = "white"))

save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/'
p_manhattan <- readRDS(paste0(save_dir,"median_CN_distance_6series.rds"))
pDLP_tree <- cowplot::plot_grid(pDLP_tree1, p_manhattan, nrow=1, rel_widths = c(1,0.2)) +
  theme(plot.background = element_rect(fill = "white", colour = "white"))

ptotal_DLP <- cowplot::plot_grid(pDLP_tree, pRxDLP, nrow=2,rel_heights = c(2.2,2))
# ptotal_DLP

save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/'
base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/'
p1 <- readRDS(paste0(save_dir,"clonealign_correlation-SA609.rds"))
res <- readRDS(paste0(save_dir,"clonealign_correlation-SA609_legends.rds"))
p1 <- p1 + labs(title='Pt4 - clonal proportions') + theme(axis.title = element_text(size=9, face="bold",family=my_font, hjust = 0.5))
# pSA609_clonealign <- cowplot::plot_grid(p1, res$colplt,res$sizeplt, ncol = 1, rel_heights = c(10,0.95,0.9))

base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/'
save_dir <- base_dir




# Clonealign Correlation for all patients
save_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/"
# stat <- data.table::fread(paste0(save_dir, '10x_dlp_stat_correlation.csv')) %>% as.data.frame()
# stat$Pt <- stat$PDX
# my_font <- "Helvetica"
# p_corr_all <- ggplot(stat, aes(x = DLP, y = clonealign)) + 
#   geom_point(aes(color = Pt, shape=Pt), alpha=0.9, size=3)  #, size=2*log10(pct_genes)size=4.5
# 
# p_corr_all <- p_corr_all + thesis_theme
# lg_pos <- "top"
# # lg_pos <- c(0.7, 0.2)
# p_corr_all <- p_corr_all + ggplot2::theme(legend.position = lg_pos,
#                         legend.title = element_blank(), 
#                         legend.key=element_blank(),
#                         legend.text = element_text(color="black",size=9, hjust = 0.5, family=my_font)) +
#   labs(x='Clonal proportion DLP', y='Clonal proportion clonealign 10x', 
#        title = '6 patients - Clonal proportions')
# 
# p_corr_all <- p_corr_all + guides(color = guide_legend(nrow = 2, override.aes = list(size=4)))
# saveRDS(p_corr_all, paste0(save_dir,'clonealign_correlation_all_patients.rds'))
save_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/"
p_corr_all <- readRDS(paste0(save_dir,'clonealign_correlation_all_patients.rds'))
# stat_DLP <- cowplot::plot_grid(res$colplt,res$colplt, ncol = 1, rel_heights = c(2, 1))
pSA609_clonealign <- cowplot::plot_grid(p_corr_all,  res$colplt, res$sizeplt + theme(legend.title = element_blank())
                                        , ncol = 1, rel_heights = c(4,0.4,0.4))
fig2_part2 <- cowplot::plot_grid(pSA609_clonealign, p1, ncol = 2, rel_widths = c(1,1))+
  theme(plot.background = element_rect(fill = "white", colour = "white"))


p_fig2 <- cowplot::plot_grid(ptotal_DLP, fig2_part2, ncol = 1, rel_heights = c(1.2,0.85), labels = c('a','b'))
ggsave(paste0(save_dir,"Fig2_DLPtree_clonealign_correlation.pdf"),
       plot = p_fig2,
       height = 11,
       width = 7.5,
       useDingbats=F,
       dpi = 150
)

ggsave(paste0(save_dir,"Fig2_DLPtree_clonealign_correlation.png"),
       plot = p_fig2,
       height = 11,
       width = 7.5,
       type = "cairo-png",
       dpi=150
)

