library(dplyr)

# SA609
script_dir <- '/home/htran/Projects/farhia_project/rnaseq/trajectory_analysis/'
source(paste0(script_dir, "slingshot_utils.R"))
source(paste0(script_dir, "tradeseq_utils.R"))
library(grid)
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/'
results_dir <- paste0(output_dir,'figs_v4/')
save_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/figs_v5/"
my_font <- "Helvetica"
datatag <- 'SA609'

## Trajectory plot
## See output of res <- plot_all_lingeages(sce, crv_umap_embed, output_dir, datatag)
l1 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage1_SA609.rds'))
l1 <- l1 + theme(plot.title = element_text(color="black", size=11, hjust = 0.5, family=my_font),
                 plot.margin = unit(c(0, 0, 0, 0), "null"),
                 panel.spacing = unit(c(0, 0, 0, 0), "null"))

l2 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage2_SA609.rds'))
l2 <- l2 + theme(plot.title = element_text(color="black", size=11, hjust = 0.5, family=my_font),
                 plot.margin = unit(c(0, 0, 0, 0), "null"),
                 panel.spacing = unit(c(0, 0, 0, 0), "null"))

l3 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage3_SA609.rds'))
l3 <- l3 + theme(plot.title = element_text(color="black", size=11, hjust = 0.5, family=my_font),
                 plot.margin = unit(c(0, 0, 0, 0), "null"),
                 panel.spacing = unit(c(0, 0, 0, 0), "null"))

total <- readRDS(paste0(results_dir,'ts_slingshot_out_wholedataset_SA609.rds'))
total <- total + theme(plot.title = element_text(color="black", size=16, hjust = 0, family=my_font, face="bold"),
                       axis.title = element_text(color="black", size=10, hjust = 0, family=my_font))
ts_color <- readRDS(paste0(results_dir,'treatment_plt.rds'))
ts_color


p1 <- cowplot::plot_grid(l3, l2, l1, ncol=1)

# p2 <- cowplot::plot_grid(total, ts_color, ncol=1, rel_heights = c(4,1))+
#   theme(plot.background = element_rect(fill = "white", colour = "white"))
p2 <- cowplot::plot_grid(total, NULL, ncol=1, rel_heights = c(3.7,1))+
  theme(plot.background = element_rect(fill = "white", colour = "white"))

ptotal_traj_609 <- cowplot::plot_grid(p2, p1, rel_widths = c(2,1), nrow = 1)

ggsave(paste0(save_dir,"hm_gene_modules_",datatag, ".png"),
       plot = ptotal_traj_609,
       height = 4,
       width = 5,
       # useDingbats=F,
       dpi=150)
ggsave(paste0(save_dir,"hm_gene_modules_",datatag, "_lg.png"),
       plot = ts_color,
       height = 3.5,
       width = 3.5,
       # useDingbats=F,
       dpi=150)


library(circlize)
library(ComplexHeatmap)
# phm_total <- cowplot::plot_grid(hm_plg, NULL, ncol=2, rel_widths = c(1,0.1),labels = c('c','d'))
# phm_total
# output_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/"
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq_v4/'
gap <- readRDS(paste0(output_dir,'Activated_Repressed_Transient_genes.rds'))
# gap

# dev.off()
gap1 <- grid.grabExpr(ComplexHeatmap::draw(gap, annotation_legend_side = "bottom",
                                           heatmap_legend_side = "bottom",
                                           padding = unit(c(1, 1, 1, 1), "mm")))

hm_plg_cistrans <- readRDS(paste0(output_dir,"gene_type_hm_legend.rds"))



save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
pw_609 <- readRDS(paste0(save_fig_dir, datatag, "_trajectory_pathways.rds"))

# lg1 <- cowplot::plot_grid(plg, hm_plg, rel_heights = c(3,1), nrow=2)
# lg1 <- cowplot::plot_grid(hm_plg_exp, hm_plg_cistrans, rel_heights = c(1,1), nrow=2)+
#   theme(plot.background = element_rect(fill = "white", colour = "white"))
# hm_lg <- cowplot::plot_grid(hm_plg_exp, hm_plg_cistrans, nrow=1)+
#   theme(plot.background = element_rect(fill = "white", colour = "white"))

hm_lg <- hm_plg_exp

hm_total <- cowplot::plot_grid(NULL, NULL, gap1, nrow=3, rel_heights = c(0.2, 0.5, 10),
                               hjust = -0.1,label_fontface="plain") + #labels=c(' ','Heatmap x: Gene Modules/ y: Pseudotime','')
  theme(plot.background = element_rect(fill = "white", colour = "white"))
# phm_total <- cowplot::plot_grid(hm_total, pgenes_ls, ncol=2, rel_widths = c(1,1),labels = c('c','d'))
ggsave(paste0(save_dir,"hm_pseudo_gene_modules_",datatag, ".png"),
       plot = hm_total,
       height = 4.5,
       width = 5,
       # useDingbats=F,
       dpi=150)

## Color plot
# ggsave(paste0(save_dir,"pseudo_gene_modules_",datatag, "_lg.png"),
#        plot = gc,
#        height = 2.5,
#        width = 3,
#        # useDingbats=F,
#        dpi=150)

ggsave(paste0(save_dir,"hm_pseudo_gene_modules_",datatag, "_lg.png"),
       plot = cowplot::ggdraw(hm_lg)+
         theme(plot.background = element_rect(fill = "white", colour = "white")),
       height = 0.5,
       width = 2,
       # useDingbats=F,
       dpi=150)

# hm_lg

# ls_pw_SA609 <- cowplot::plot_grid(pw_609, hm_lg, rel_widths = c(3,1), nrow=1)+ 
#   theme(plot.background = element_rect(fill = "white", colour = "white"))
# phm_total <- cowplot::plot_grid(hm_total, ls_pw_SA609, ncol=1, rel_heights = c(3,1),labels = c('b','c'))

# p609 <- cowplot::plot_grid(ptotal, phm_total, nrow=2, rel_heights = c(1,1.7), labels=c('a',''))
p609 <- cowplot::plot_grid(ptotal_traj_609, hm_total, nrow=2, rel_heights = c(1,1.5), labels=c('a','b'))
png(paste0(results_dir,"Fig7_trajectory_",datatag,".png"), height = 2*1200, width=2*950,res = 2*72)
print(p609)
dev.off()



## SA535
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/'
results_dir <- paste0(output_dir,'figs_v4/')
my_font <- "Helvetica"
datatag <- 'SA535'

# library(extrafont)
# font_import(prompt=F, paths ='/usr/share/fonts/truetype/myfonts/') # import Helvetica font
# fonts()

l1 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage1_',datatag,'.rds'))
l1 <- l1 + theme(plot.title = element_text(color="black", size=11, hjust = 0.5, family=my_font),
plot.margin = unit(c(0, 0, 0, 0), "null"),
panel.spacing = unit(c(0, 0, 0, 0), "null"))

l2 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage2_',datatag,'.rds'))
l2 <- l2 + theme(plot.title = element_text(color="black", size=11, hjust = 0.5, family=my_font),
plot.margin = unit(c(0, 0, 0, 0), "null"),
panel.spacing = unit(c(0, 0, 0, 0), "null"))

l3 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage3_',datatag,'.rds'))
l3 <- l3 + theme(plot.title = element_text(color="black", size=11, hjust = 0.5, family=my_font),
plot.margin = unit(c(0, 0, 0, 0), "null"),
panel.spacing = unit(c(0, 0, 0, 0), "null"))

l4 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage4_',datatag,'.rds'))
l4 <- l4 + theme(plot.title = element_text(color="black", size=11, hjust = 0.5, family=my_font),
plot.margin = unit(c(0, 0, 0, 0), "null"),
panel.spacing = unit(c(0, 0, 0, 0), "null"))

total <- readRDS(paste0(results_dir,'ts_slingshot_out_wholedataset_',datatag,'.rds'))
total <- total + theme(plot.title = element_text(color="black", size=16, hjust = 0, family=my_font, face="bold"),
                       axis.title = element_text(color="black", size=10, hjust = 0.5, family=my_font))

ts_color <- readRDS(paste0(results_dir,'treatment_plt.rds'))
# ts_color <- pc
# lg <- pc
p1 <- cowplot::plot_grid(l4,l2, ncol=1) + theme(plot.background = element_rect(fill = "white", colour = "white"))
# p2 <- cowplot::plot_grid(l1, l3, ts_color, ncol=3)+ 
#   theme(plot.background = element_rect(fill = "white", colour = "white"))

p2 <- cowplot::plot_grid(NULL, l3, l1, ncol=3)+ 
  theme(plot.background = element_rect(fill = "white", colour = "white"))
p3 <- cowplot::plot_grid(total, p1, rel_widths = c(2,1), nrow = 1)
ptotal_traj_SA535 <- cowplot::plot_grid(p3, p2, rel_heights = c(2,1), ncol=1)+ 
  theme(plot.background = element_rect(fill = "white", colour = "white"))
ptotal_traj_SA535
dev.off()


ggsave(paste0(save_dir,"pseudo_gene_modules_",datatag, ".png"),
       plot = ptotal_traj_SA535,
       height = 4,
       width = 5,
       # useDingbats=F,
       dpi=150)

ggsave(paste0(save_dir,"pseudo_gene_modules_",datatag, "_color_lg.png"),
       plot = ts_color,
       height = 2.5,
       width = 3.5,
       # useDingbats=F,
       dpi=150)

# Genes module heatmap
gap_SA535 <- readRDS(paste0(output_dir,'Activated_Repressed_Transient_genes.rds'))
# gap_SA535

hm_plg_cistrans <- readRDS(paste0(output_dir,"gene_type_hm_legend.rds"))
col_fun = colorRamp2(c(-4,0,4), c("blue", "white", "red"))
lgd = ComplexHeatmap::Legend(col_fun = col_fun, title = "Avg Exp", at = c(-4,-2,0, 2, 4), direction = "horizontal")
hm_plg_exp <- grid.grabExpr(ComplexHeatmap::draw(lgd))
lg1 <- cowplot::plot_grid(hm_plg_exp, hm_plg_cistrans, rel_heights = c(1,1), nrow=2)
lg1 <- hm_plg_cistrans

save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
pw_535 <- readRDS(paste0(save_fig_dir, datatag, "_trajectory_pathways.rds"))
# pw_535


gap1_SA535 <- grid.grabExpr(ComplexHeatmap::draw(gap_SA535, annotation_legend_side = "bottom",
                                           heatmap_legend_side = "bottom",
                                           padding = unit(c(5, 1, 1, 1), "mm")))

hm_total_SA535 <- cowplot::plot_grid(NULL, NULL, gap1_SA535, nrow=3, rel_heights = c(0.2, 0.5, 10),
                               hjust = -0.1,label_fontface="plain") + #labels=c(' ','Heatmap x: Gene Modules/ y: Pseudotime','')
  theme(plot.background = element_rect(fill = "white", colour = "white"))

ggsave(paste0(save_dir,"hm_pseudo_gene_modules_",datatag, "_hm.png"),
       plot = hm_total_SA535,
       height = 4.8,
       width = 5,
       # useDingbats=F,
       dpi=150)

ggsave(paste0(save_dir,"hm_pseudo_gene_modules_",datatag, "_lg.png"),
       plot = cowplot::ggdraw(hm_plg_cistrans)+
         theme(plot.background = element_rect(fill = "white", colour = "white")),
       height = 1.5,
       width = 1,
       # useDingbats=F,
       dpi=150)

ls_pw_SA535 <- cowplot::plot_grid(pw_535, hm_plg_cistrans, rel_widths = c(3,1), nrow=1)+ 
  theme(plot.background = element_rect(fill = "white", colour = "white"))
# phm_total_SA535 <- cowplot::plot_grid(hm_total_SA535, ls_pw_SA535, ncol=1, rel_heights = c(3,1),labels = c('e','f'))



# p535 <- cowplot::plot_grid(ptotal_traj_SA535, phm_total_SA535, nrow=2, rel_heights = c(1,1.7), labels = c('d',''))
p535 <- cowplot::plot_grid(ptotal_traj_SA535, hm_total_SA535, nrow=2, rel_heights = c(1,1.5), labels = c('c','d'))

final_plt <- cowplot::plot_grid(p609, p535, ncol=2, rel_widths = c(1,1))

final_plt2 <- cowplot::plot_grid(final_plt, pathway_all, ncol=1, rel_heights = c(2,0.9))
png(paste0(results_dir,"Fig7_trajectory_SA609_SA535.png"), height = 2*1350, width=2*1000,res = 2*72)
print(final_plt2)
dev.off()


ggsave(paste0(results_dir,"Fig7_trajectory_SA609_SA535.png"),
       plot = final_plt2,
       height = 16,
       width = 13,
       # useDingbats=F,
       dpi=150)

ggsave(paste0(results_dir,"Fig7_trajectory_SA609_SA535.pdf"),
       plot = final_plt2,
       height = 16,
       width = 13,
       useDingbats=F,
       dpi=150)
