

output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/'
results_dir <- paste0(output_dir,'figs_v3/')
my_font <- "Helvetica"
datatag <- 'SA1035'
library(extrafont)
font_import(prompt=F, paths ='/usr/share/fonts/truetype/myfonts/') # import Helvetica font
fonts()

l1 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage1_',datatag,'.rds'))
# l1 <- l1 + theme(plot.title = element_text(color="black", size=13, hjust = 0.5, family=my_font))

l2 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage2_',datatag,'.rds'))
# l2 <- l2 + theme(plot.title = element_text(color="black", size=12, hjust = 0.5, family=my_font))

l3 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage3_',datatag,'.rds'))
# l3 <- l3 + theme(plot.title = element_text(color="black", size=12, hjust = 0.5, family=my_font))

l4 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage4_',datatag,'.rds'))
# l4 <- l4 + theme(plot.title = element_text(color="black", size=12, hjust = 0.5, family=my_font))

l5 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage4_',datatag,'.rds'))

total <- readRDS(paste0(results_dir,'ts_slingshot_out_wholedataset_',datatag,'.rds'))
total <- total + theme(plot.title = element_text(color="black", size=20, hjust = 0, family=my_font, face="bold"))
ts_color <- readRDS(paste0(results_dir,'treatment_plt.rds'))

dev.off()
p1 <- cowplot::plot_grid(l3,l2, ncol=1)
p2 <- cowplot::plot_grid(l4, l1, ts_color, ncol=3) +
  theme(plot.background = element_rect(fill = "white", colour = "white"))
p3 <- cowplot::plot_grid(total, p1, rel_widths = c(2,1), nrow = 1)
ptotal <- cowplot::plot_grid(p3, p2, rel_heights = c(2,1), ncol=1)
ptotal

png(paste0(results_dir,"slingshot_out_summary_",datatag,"_part11.png"), height = 2*500, width=2*700,res = 2*72)
print(ptotal)
dev.off()


pts <- readRDS(paste0(results_dir,datatag,'_treatment_desc_trajectory_prevalence.rds'))
pc <- readRDS(paste0(results_dir,datatag,'_treatment_clone_trajectory_prevalence.rds'))
# gb_heatmap = grid.grabExpr(draw(Heatmap(mat)))
prevalence_ts <- grid.grabExpr(ComplexHeatmap::draw(pts, annotation_legend_side = "bottom",
                                                    heatmap_legend_side = "bottom",
                                                    padding = unit(c(8, 6, 4, 2), "mm")))
prevalence_clone <- grid.grabExpr(ComplexHeatmap::draw(pc, annotation_legend_side = "bottom",
                                                       heatmap_legend_side = "bottom",
                                                       padding = unit(c(8, 2, 4, 2), "mm")))
pvalenc <- cowplot::plot_grid(prevalence_ts, prevalence_clone,rel_widths = c(1,1.4), ncol=2)
ptraj <- cowplot::plot_grid(ptotal, pvalenc,
                            rel_widths = c(3,2), ncol=2, labels = c('a','b'))
png(paste0(results_dir,"trajectory_",datatag,"_part1.png"), height = 2*400, width=2*950,res = 2*72)
print(ptraj)
dev.off()






datatag <- 'SA1035'
results_dir <- paste0("/home/htran/storage/datasets/drug_resistance/rna_results/",datatag,"_rna/slingshot_trajectory/figs_v3/")
demo_genes <- c('S100A6','CST3','NFKBIA','CSTB')
titles <- c('S100A6: M5','CST3: M3','NFKBIA: M4','CSTB: M2')

plt_genes_demo <- list()

## First gene: plot legend
i = 1
g <- demo_genes[i]
exg <- readRDS(paste0(results_dir,g,'.rds'))  
exg <- exg + theme(plot.title = element_text(color="black", size=13, hjust = 0.5, family=my_font),
                   legend.text = element_text(color="black", size=10, hjust = 0.5, family=my_font))
exg <- exg + labs(title = titles[i])
plt_genes_demo[[g]] <- exg

i = 2
g <- demo_genes[i]
exg <- readRDS(paste0(results_dir,g,'_without_lg.rds'))  
exg <- exg + theme(plot.title = element_text(color="black", size=13, hjust = 0.5, family=my_font))
exg <- exg + labs(title = titles[i])
plt_genes_demo[[g]] <- exg

i = 3
g <- demo_genes[i]
exg <- readRDS(paste0(results_dir,g,'_without_lg.rds'))  
exg <- exg + theme(plot.title = element_text(color="black", size=13, hjust = 0.5, family=my_font))
exg <- exg + labs(title = titles[i])
plt_genes_demo[[g]] <- exg

i = 4
g <- demo_genes[i]
exg <- readRDS(paste0(results_dir,g,'_without_lg.rds'))  
exg <- exg + theme(plot.title = element_text(color="black", size=13, hjust = 0.5, family=my_font))
exg <- exg + labs(title = titles[i])
plt_genes_demo[[g]] <- exg

dev.off()
pgene <- cowplot::plot_grid(plotlist = plt_genes_demo, ncol=2) #, labels = c('a','b','c','d')
pgene


library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(c(-4,0,4), c("blue", "white", "red"))
lgd = ComplexHeatmap::Legend(col_fun = col_fun, title = "Avg Exp", at = c(-4,-2,0, 2, 4), direction = "horizontal")
hm_plg_exp <- grid.grabExpr(ComplexHeatmap::draw(lgd))
hm_plg_cistrans <- readRDS(paste0(output_dir,"gene_type_hm_legend.rds"))
lg1 <- cowplot::plot_grid(hm_plg_exp, hm_plg_cistrans, rel_heights = c(1,1), nrow=2)+
  theme(plot.background = element_rect(fill = "white", colour = "white"))

gap <- readRDS(paste0(dirname(results_dir),'/Activated_Repressed_Transient_genes.rds'))
gap1 <- grid.grabExpr(ComplexHeatmap::draw(gap, annotation_legend_side = "bottom",
                                           heatmap_legend_side = "bottom",
                                           padding = unit(c(11, 2, 3, 1), "mm")))





## Pathway analysis
# save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
# pathway_stat <- data.table::fread(paste0(save_fig_dir,'pw_stats_modules_3series.csv')) %>% as.data.frame()
# pathway_stat <- pathway_stat %>%
#   dplyr::filter(datatag=='SA1035') %>%
#   dplyr::rename(pathway=reference_set, nb_signf_genes=nb_signif_genes)
# pathway_stat$gene_type_module <- gsub('Module','M',pathway_stat$gene_type_module)
# my_font <- "Helvetica"
# pathway_stat$pathway <- tolower(gsub('HALLMARK_','',pathway_stat$pathway))
# pathway_stat$datatag <- pathway_stat$gene_type_module
# pathway_stat$pathway <- gsub('_response','',pathway_stat$pathway)
# pw_1035 <- viz_pathways(pathway_stat, datatag, save_fig_dir)
# pw_1035 <- readRDS(paste0(save_fig_dir, datatag, "_trajectory_pathways.rds"))

save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
pw_1035 <- readRDS(paste0(save_fig_dir, datatag, "_trajectory_pathways.rds"))
pw_1035




pgene_lg <- cowplot::plot_grid(lg1, pw_1035, rel_widths = c(1,3.8), nrow=1)
pgenes_ls <- cowplot::plot_grid(pgene, pgene_lg, ncol=1, rel_heights = c(2,1))

hm_total <- cowplot::plot_grid(NULL, NULL, gap1, nrow=3, rel_heights = c(0.2, 0.5, 10),
                               hjust = -0.1,label_fontface="plain",
                               labels=c(' ','Heatmap  x: Gene Modules/ y: Pseudotime','')) + 
  theme(plot.background = element_rect(fill = "white", colour = "white"))

phm_total <- cowplot::plot_grid(hm_total, pgenes_ls, ncol=2, rel_widths = c(1,1),labels = c('c','d'))

p1035 <- cowplot::plot_grid(ptraj, phm_total, nrow=2, rel_heights = c(1,1.5))


ggsave(paste0(results_dir,"SUPP_Fig_trajectory_",datatag,".png"),
       plot = p1035,
       height = 13.5,
       width = 12.2,
       # useDingbats=F,
       dpi=150)


ggsave(paste0(results_dir,"SUPP_Fig_trajectory_",datatag,".pdf"),
       plot = p1035,
       height = 13.5,
       width = 12.2,
       useDingbats=F,
       dpi=150)


