# Figure 2: DLP tree, clonealign 


# base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/'
# color_df <- data.table::fread(paste0(base_dir,'figs/colorcode_total.csv'))
# dim(color_df)
# head(color_df)
# color_df <- color_df %>%
#   dplyr::filter(datatag!='SA1035')
# 
# cols_use <- c('H'='#521913','F'='#ABB9ED','G'='#C1221F','I'='#E5702F',
#              'K'='#A1BE58','E'='#4F97BA','C'='#65AD97','A'='#5666B6','D'='#92559E',
#              'J'='#BFA0CC','B'='#DDA83B')
# 
# length(cols_use)
# color_df1 <- data.frame(clone_id=names(cols_use),colour=cols_use)
# color_df1$datatag <- 'SA1035'
# color_df <- dplyr::bind_rows(color_df, color_df1)
# summary(as.factor(color_df$datatag))
# data.table::fwrite(color_df, paste0(base_dir,'figs/colorcode_total_v2.csv'))

library(dplyr)
library(ggplot2)
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
                            labels = c('Pt1','Pt2','Pt3'), label_y = 0.2, hjust = -0.5)
pSeries <- cowplot::plot_grid(dlp_609,dlp_535, dlp_1035, ncol=3, rel_widths = c(0.8,1,1),
                              labels = c('Pt4','Pt5','Pt6'), label_y = 0.2, hjust = -0.5)
ptotal_DLP <- cowplot::plot_grid(NULL, pUnRx,NULL, pSeries, nrow=4, rel_heights = c(0.1,1,0.1,1),
                             labels = c('     -Rx PDX tumors','',
                                        '     -Rx/+Rx time series PDX tumors',''),
                             # , label_x = 0, label_y = 1
                             label_x = .01, hjust = 0, label_size=13)
ptotal_DLP
save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/'
png(paste0(save_dir,"Fig2_v2_part1.png"), height = 2*350, width=2*650, res = 2*60)
print(ptotal_DLP)
dev.off()

ggsave(paste0(save_dir,"Fig2_v2_part1.png"),
       plot = ptotal_DLP,
       height = 3.8,
       width = 6.5,
       type = "cairo-png",
       dpi=150
)


base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/'
p1 <- readRDS(paste0(save_dir,"clonealign_correlation-SA609.rds"))
res <- readRDS(paste0(save_dir,"clonealign_correlation-SA609_legends.rds"))
p1 <- p1 + labs(title='Pt4 - clonal proportions')
pSA609_clonealign <- cowplot::plot_grid(p1, res$colplt,res$sizeplt, ncol = 1, rel_heights = c(10,0.95,0.9))

save_dir <- base_dir
p_corr <- readRDS(paste0(save_dir,"clonealign_correlation-all_v2.rds"))
ggsave(paste0(save_dir,"Fig2_v2_part2.png"),
       plot = p_corr,
       height = 4.2,
       width = 6.5,
       type = "cairo-png",
       dpi=150
)

pfig2 <- cowplot::plot_grid(ptotal_DLP, p_corr, ncol = 1, rel_heights = c(1,1.3), labels = c('a','b'))

png(paste0(save_dir,"Fig2_DLPtree_clonealign_correlation.png"), height = 2*800, width=2*650, res = 2*72)
print(pfig2)
dev.off()

ggsave(paste0(save_dir,"Fig2_DLPtree_clonealign_correlation.pdf"),
              plot = pfig2,
              height = 8.5,
              width = 7.5,
              useDingbats=F,
              dpi = 300
       )

ggsave(paste0(save_dir,"Fig2_v2.png"),
       plot = pfig2,
       height = 8.5,
       width = 8,
       type = "cairo-png",
       dpi=200
)
