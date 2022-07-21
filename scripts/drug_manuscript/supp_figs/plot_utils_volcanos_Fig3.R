
library(dplyr)
library(cowplot)


# Utils function here
source('/home/htran/Projects/farhia_project/rnaseq/pipeline/utils/plot_utils.R')
source('/home/htran/Projects/farhia_project/rnaseq/cis_trans/in_cis_trans_utils.R')
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'

# output_dir <- paste0(base_dir,'SA609_rna/deg_analysis/SA609-v6/SA609_UTTTT_R_UUUUU_H/')
# de_genes_609 <- read.csv(paste0(output_dir,'signif_genes.csv'), check.names = F, stringsAsFactors = F)
# 
# summary(as.factor(de_genes_609$Gene_Type))
input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/signif_genes/'
deg_fn <- paste0(input_dir,'scrande_SA609_4/signif_genes.csv')
de_genes_609 <- data.table::fread(deg_fn) %>% as.data.frame()
dim(de_genes_609)
datatag <- 'SA609'


# plttitle <- 'SA609: Rx X7:cloneA vs. UnRx X7:cloneH'
plttitle <- 'Pt4: Rx X7:cloneA vs. UnRx X7:cloneH'
# de_genes_609 <- de_genes_609 %>%
#   dplyr::filter(abs(logFC)>0.5)
# vol_609 <- plot_DE_genes_edgeR(de_genes, topGenes, capstr='', FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
#                     plttitle, save_dir=output_dir, legendVisible=F,
#                     iscaption=TRUE, legend_verbose='none', save_plot=TRUE)
# vol_609

vol_609 <- plot_DE_genes_ggplot_Fig3(de_genes_609, NULL, capstr='', 
                                FDRcutoff=0.01, logFCcutoff=0.5, pValuecutoff=0.05,
                                plttitle = NULL, save_dir=output_dir, legendVisible=F,
                                iscaption=TRUE, legend_verbose='none', save_plot=F,
                                xl=c(-3.5,3.5), yl=c(0, 320))

vol_609_legend <- plot_DE_genes_ggplot_Fig3(de_genes_609, NULL, 
                                       capstr='', FDRcutoff=0.01, logFCcutoff=0.5, pValuecutoff=0.05,
                                       plttitle = NULL, save_dir=output_dir, legendVisible=F,
                                       iscaption=TRUE, legend_verbose='left', save_plot=F)
saveRDS(vol_609, paste0(save_fig_dir,'vol_609.rds'))
saveRDS(vol_609_legend, paste0(save_fig_dir,'vol_609_legend.rds'))

vol_609$p_prop_trans
vol_609$p_vol_cis
vol_609$p_vol_trans
vol_609$p_prop_cis
vol_609$p_prop_trans
volplotls <- cowplot::plot_grid(vol_609$p_vol_trans,vol_609$p_vol_cis,
                                ncol=2, align = 'v') + theme(plot.background = element_rect(fill = "white", colour = "white"))
volplot_prop <- cowplot::plot_grid(vol_609$p_prop_trans, vol_609_legend$plg_trans, 
                                   vol_609$p_prop_cis, vol_609_legend$plg_cis,
                                    nrow=1, align = 'h', rel_widths = c(0.5,2,0.5,2)
                                    )+
                                    theme(plot.background = element_rect(fill = "white", colour = "white"))
volplot_prop
# volplot_prop <- cowplot::plot_grid(vol_609_legend$plg_trans, vol_609_legend$plg_cis,
#                                    nrow=1, align = 'v',
#                                    rel_widths = c(1,1))+
#   theme(plot.background = element_rect(fill = "white", colour = "white"))
# ptotal <- cowplot::plot_grid(volplotls, volplot_prop, rel_widths = c(2.5,1.5))

# ,align = 'v'
# volplotls <- cowplot::plot_grid(vol_609$p_vol_trans,vol_609$p_vol_cis, 
#                                 vol_609_legend$plg_trans, vol_609_legend$plg_cis,
#                                 ncol=2, rel_heights = c(2,2,1,1)) + theme(plot.background = element_rect(fill = "white", colour = "white"))
# volplot_prop <- cowplot::plot_grid(vol_609$p_prop_trans, vol_609_legend$plg_trans, NULL,
#                                    vol_609$p_var,NULL,vol_609$p_prop_cis, vol_609_legend$plg_cis, 
#                                    ncol=1, align = 'v', 
#                                    rel_heights = c(1,1.5,0.3,2.5,0.3,1,2.5))+
#   theme(plot.background = element_rect(fill = "white", colour = "white"))
# volplot_prop
# ptotal <- cowplot::plot_grid(volplotls, volplot_prop, rel_widths = c(2.5,1.5))
ptotal <- cowplot::plot_grid(volplotls, volplot_prop, rel_heights = c(2,1), ncol=1)
save_fig_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/trackplots/"
track_609 <- readRDS(paste0(save_fig_dir,'main_plot_SA609.rds'))
# track_609 <- main_plot_SA609
# df <- mtcars[, c("mpg", "cyl", "wt")]
# df$cyl <- as.factor(df$cyl)
# pt <- ggplot(df, aes(x=wt, y=mpg)) +
#   geom_point(shape=18)
ptotal_vol <- cowplot::plot_grid(ptotal, pdis_logFC,
                                       ncol=1, rel_heights = c(1,1), labels = c('b','c')
                                       ) +#
  theme(plot.background = element_rect(fill = "white", colour = "white"))

# png(paste0(save_fig_dir,"Fig4_volcanos_track_SA609.png"), height = 2*900, width=2*650, res = 2*72)
# print(ptotal_track_vol)
# dev.off()

ggsave(paste0(save_fig_dir,"Fig3_volcanos_track_SA609.pdf"),
       plot = ptotal_track_vol,
       height = 14,
       width = 10,
       useDingbats=F,
       dpi = 300)

ggsave(paste0(save_fig_dir,"Fig3_volcanos_track_SA609.png"),
       plot = ptotal_track_vol,
       height = 14,
       width = 10,
       # useDingbats=F,
       type = "cairo-png",
       dpi = 300
)

ggsave(paste0(save_fig_dir,"Fig3_volcanos_track_SA609_p1.png"),
       plot = track_609,
       height = 6,
       width = 10,
       # useDingbats=F,
       type = "cairo-png",
       dpi = 300
)

ggsave(paste0(save_fig_dir,"Fig3_volcanos_track_SA609_p2.png"),
       plot = ptotal_vol,
       height = 8,
       width = 10,
       # useDingbats=F,
       type = "cairo-png",
       dpi = 300
)


# vol_609 <- plot_DE_genes_edgeR(de_genes, topGenes, capstr='', 
#                                FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
#                                plttitle, output_dir, legendVisible=F,
#                                iscaption=TRUE)

# cis_6 <- de_genes %>%
#   dplyr::filter(!is.na(classified_gene_dlp))
# 
# trans_6 <- de_genes %>%
#   dplyr::filter(is.na(classified_gene_dlp))
# 
# round(mean(abs(cis_6$logFC)),2) #0.62
# round(sd(abs(cis_6$logFC)),2)   #0.35
# round(mean(abs(trans_6$logFC)),2) #0.53
# round(sd(abs(trans_6$logFC)),2) #0.29
# stat_609 <- compute_stat(de_genes_609, 'SA609')
# stat_609_nb <- compute_stat_nb(de_genes_609, 'SA609')
# stat_609_top <- compute_stat_top50(de_genes_609, 'SA609')
# cis_609_top <- compute_stat_top50(de_genes_609, 'in-cis')
# top5pct_609 <- compute_stat_top5percent(de_genes_609, 'SA609')

# dlp_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/clonealign/')
# save_dir <- paste0(dirname(dlp_dir),'/cis_trans/')
# output_dir <- paste0(save_dir,'rstmm6_SA1035_UTTTT_H_UUUUU_E/')
# de_genes_1035 <- data.table::fread(paste0(output_dir,'signif_genes_FDR0.01.csv')) %>% as.data.frame()
datatag <- 'SA1035'
tag <- 'Pt6'
# Utils function here 
source('/home/htran/Projects/farhia_project/rnaseq/cis_trans/in_cis_trans_utils.R')
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'

de_genes_1035 <- data.table::fread(paste0(output_dir,datatag,'/SA1035_2_SA1035_UTTTT_H_UUUUU_E/signif_genes_FDR_0.01.csv')) %>% as.data.frame()
dim(de_genes_1035)
# plttitle <- 'SA1035: Rx X8:cloneH vs UnRx X8:cloneE'
plttitle <- paste0(tag, ': Rx X8:cloneH vs UnRx X8:cloneE')

save_fig_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/trackplots/"
track_1035 <- readRDS(paste0(save_fig_dir,'main_plot_SA1035.rds'))

vol_1035 <- plot_DE_genes_ggplot(de_genes_1035, get_top_genes(de_genes_1035), capstr='', 
                                 FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                                 plttitle, save_dir=output_dir, legendVisible=T,
                                 iscaption=TRUE, legend_verbose='none', save_plot=F,
                                 xl=c(-3.5,3.5))
vol_1035$p_vol_cis
vol_1035$p_vol_trans
vol_1035_legend <- plot_DE_genes_ggplot(de_genes_1035, get_top_genes(de_genes_1035), 
                                        capstr='', FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                                        plttitle, save_dir=output_dir, legendVisible=F,
                                       iscaption=TRUE, legend_verbose='top', save_plot=F)
saveRDS(vol_1035, paste0(save_fig_dir,'vol_',datatag,'.rds'))
saveRDS(vol_1035_legend, paste0(save_fig_dir,'vol_',datatag,'_legend.rds'))
plg_cis <- cowplot::ggdraw() + cowplot::draw_plot(cowplot::get_legend(vol_1035_legend$p_vol_cis +  guides(color = guide_legend(nrow = 2, title = NULL))))
plg_trans <- cowplot::ggdraw() + cowplot::draw_plot(cowplot::get_legend(vol_1035_legend$p_vol_trans +  guides(color = guide_legend(nrow = 1, title = NULL))))

volplotls <- cowplot::plot_grid(vol_1035$p_vol_cis, vol_1035$p_vol_trans,
                                ncol=1, align = 'v')
volplot_prop <- cowplot::plot_grid(vol_1035$p_prop_cis, plg_cis, 
                                   vol_1035$p_prop_trans, plg_trans, 
                                   ncol=1, align = 'v', 
                                   rel_heights = c(3,1.5,3,1))

ptotal <- cowplot::plot_grid(volplotls, volplot_prop, rel_widths = c(3,2))
ptotal_track_vol <- cowplot::plot_grid(track_1035, ptotal, rel_heights = c(1,1.2), ncol=1)

png(paste0(save_fig_dir,"SUPPFig_volcanos_track_",datatag,".png"), height = 2*900, width=2*650, res = 2*72)
print(ptotal_track_vol)
dev.off()

ggsave(paste0(save_fig_dir,"SUPPFig_volcanos_track_",datatag,".pdf"),
       plot = ptotal_track_vol,
       height = 9.5,
       width = 8,
       useDingbats=F,
       dpi = 300)

# 
# cis_10 <- de_genes %>%
#   dplyr::filter(!is.na(classified_gene_dlp))
# 
# trans_10 <- de_genes %>%
#   dplyr::filter(is.na(classified_gene_dlp))
# 
# round(mean(abs(cis_10$logFC)),2) #0.96
# round(sd(abs(cis_10$logFC)),2)   #0.5
# 
# round(mean(abs(trans_10$logFC)),2) #0.84
# round(sd(abs(trans_10$logFC)),2) #0.48
# 
# stat_1035 <- compute_stat(de_genes, 'SA1035')
# stat_1035_nb <- compute_stat_nb(de_genes, 'SA1035')
# stat_1035_top <- compute_stat_top50(de_genes, 'SA1035')
# cis_1035_top <- compute_stat_top50(de_genes, 'in-cis')
# top5pct_1035 <- compute_stat_top5percent(de_genes_1035, 'SA1035')


# vol_1035 <- plot_DE_genes_edgeR(de_genes, topGenes, capstr='', FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
#                     plttitle, save_dir=output_dir, legendVisible=T,
#                     iscaption=TRUE, legend_verbose='right', save_plot=TRUE)

# vol_1035 <- plot_DE_genes_ggplot(de_genes, topGenes, capstr='', FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
#                                 plttitle, save_dir=output_dir, legendVisible=T,
#                                 iscaption=TRUE, legend_verbose='none', save_plot=TRUE,xl=xl)
# vol_1035 <- plot_DE_genes_edgeR(de_genes, topGenes, capstr='', 
#                                 FDRcutoff=0.01, logFCcutoff=0.25, 
#                                 plttitle, output_dir, legendVisible=T,
#                                 iscaption=TRUE)



datatag <- 'SA535'
tag <- 'Pt5'
# dlp_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/clonealign/')
# save_dir <- paste0(dirname(dlp_dir),'/cis_trans/')
# output_dir <- paste0(save_dir,'rstmm10_SA535_UUTTTT_A_UUUUUU_G/')
# de_genes_535cis <- data.table::fread(paste0(output_dir,'signif_genes_FDR0.01.csv')) %>% as.data.frame()
save_fig_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/dlp_cnv/trackplots/"
track_535 <- readRDS(paste0(save_fig_dir,'main_plot_SA535.rds'))
de_genes_535cis <- data.table::fread(paste0(output_dir,datatag,'/SA535_2_SA535_UUTTTT_A_UUUUUU_G/signif_genes_FDR_0.01.csv')) %>% as.data.frame()

dim(de_genes_535cis)
max(de_genes_535cis$logFC)
# sum(abs(de_genes_535cis$logFC)>0.5)
# de_genes_535cis <- de_genes_535cis %>%
#   dplyr::filter(abs(logFC)>0.5)
# dim(de_genes_535cis)
# summary(as.factor(de_genes_535cis$Gene_Type))
# plttitle <- 'SA535:Cisplatin: Rx X10:cloneT vs UnRx X9:cloneJ' # combined version
# plttitle <- 'SA535:Cisplatin: Rx X10:cloneA vs UnRx X9:cloneG'
plttitle <- paste0(tag, ': Rx X10:cloneA vs UnRx X9:cloneG')
vol_535 <- plot_DE_genes_ggplot(de_genes_535cis, get_top_genes(de_genes_535cis), capstr='', 
                                 FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                                 plttitle, save_dir=output_dir, legendVisible=T,
                                 iscaption=TRUE, legend_verbose='none', save_plot=F,
                                 xl=c(-3.5,3.5))
vol_535$p_vol_cis
vol_535$p_vol_trans
vol_535_legend <- plot_DE_genes_ggplot(de_genes_535cis, get_top_genes(de_genes_535cis), 
                                        capstr='', FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                                        plttitle, save_dir=output_dir, legendVisible=F,
                                        iscaption=TRUE, legend_verbose='top', save_plot=F)
saveRDS(vol_535, paste0(save_fig_dir,'vol_',datatag,'.rds'))
saveRDS(vol_535_legend, paste0(save_fig_dir,'vol_',datatag,'_legend.rds'))
plg_cis <- cowplot::ggdraw() + cowplot::draw_plot(cowplot::get_legend(vol_535_legend$p_vol_cis +  guides(color = guide_legend(nrow = 2, title = NULL))))
plg_trans <- cowplot::ggdraw() + cowplot::draw_plot(cowplot::get_legend(vol_535_legend$p_vol_trans +  guides(color = guide_legend(nrow = 1, title = NULL))))

volplotls <- cowplot::plot_grid(vol_535$p_vol_cis, vol_535$p_vol_trans,
                                ncol=1, align = 'v')
volplot_prop <- cowplot::plot_grid(vol_535$p_prop_cis, plg_cis, 
                                   vol_535$p_prop_trans, plg_trans, 
                                   ncol=1, align = 'v', 
                                   rel_heights = c(3,1.5,3,1))

ptotal <- cowplot::plot_grid(volplotls, volplot_prop, rel_widths = c(3,2))
ptotal_track_vol <- cowplot::plot_grid(track_535, ptotal, rel_heights = c(1,1.2), ncol=1)

png(paste0(save_fig_dir,"SUPPFig_volcanos_track_",datatag,".png"), height = 2*900, width=2*650, res = 2*72)
print(ptotal_track_vol)
dev.off()

ggsave(paste0(save_fig_dir,"SUPPFig_volcanos_track_",datatag,".pdf"),
       plot = ptotal_track_vol,
       height = 9.5,
       width = 8,
       useDingbats=F,
       dpi = 300)


# meta <- data.frame(logFC_025=c(5612,3634,3973),logFC_05=c(3063,3018,2137))
# rownames(meta) <- c('SA609','SA1035','SA535')
# View(meta)

output_dir <- paste0(base_dir,'SA535_total_rna_v2/SA535-v6/SA535_UUTTTT_T_UUUUUU_J/')
# plttitle <-'SA535 cisplatin: UUTTTT-T vs UUUUUU-J'
# plttitle <-'SA535 cisplatin: X10-Rx T vs X9-UnRx J'
rm(de_genes)
de_genes_535cis <- read.csv(paste0(output_dir,'signif_genes.csv'), check.names = F, stringsAsFactors = F)
de_genes <- de_genes_535cis
dim(de_genes)
# summary(as.factor(de_genes$Gene_Type))
# cis_5c <- de_genes %>%
#   dplyr::filter(!is.na(classified_gene_dlp))
# 
# trans_5c <- de_genes %>%
#   dplyr::filter(is.na(classified_gene_dlp))
# 
# # 0.57  0.34 vs 0.54  0.31
# round(mean(abs(cis_5c$logFC)),2) #0.57
# round(sd(abs(cis_5c$logFC)),2)   #0.34
# 
# round(mean(abs(trans_5c$logFC)),2) #0.54
# round(sd(abs(trans_5c$logFC)),2) #0.31
# 
# stat_5c <- compute_stat(de_genes, 'SA535_cisplatin')
# stat_5c_nb <- compute_stat_nb(de_genes, 'SA535_cisplatin')
# 
# stat_5c_top <- compute_stat_top50(de_genes, 'SA535_cisplatin')
# cis_5c_top <- compute_stat_top50(de_genes, 'in-cis')
# top5pct_535cis <- compute_stat_top5percent(de_genes_535cis, 'SA535_cisplatin')

xl <- c(min(de_genes$logFC),3.5)
# vol_535_cis <- plot_DE_genes_edgeR(de_genes, topGenes, capstr='', FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
#                                   plttitle, save_dir=output_dir, legendVisible=F,
#                                   iscaption=TRUE, legend_verbose='none', save_plot=TRUE)
# vol_535_cis <- plot_DE_genes_edgeR(de_genes, topGenes, capstr='', 
#                                    FDRcutoff=0.01, logFCcutoff=0.25,
#                                    plttitle, output_dir, legendVisible=F,
#                                    iscaption=TRUE)
vol_535_cis <- plot_DE_genes_ggplot(de_genes, topGenes, capstr='', FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                                   plttitle, save_dir=output_dir, legendVisible=F,
                                   iscaption=TRUE, legend_verbose='none', save_plot=TRUE, xl=xl)
# ggdraw(vol_535_cis$p_vol)



output_dir <- paste0(base_dir,'SA535_total_rna_v2/SA535-v6/SA535_UXXXX_U_UUUUU_J/')
# plttitle <-'SA535 CX5461: UXXXX-U vs UUUUU-J'
# plttitle <-'SA535 CX5461: X8-Rx U vs X9-UnRx J'
plttitle <- 'SA535:CX5461: Rx X8:cloneU vs UnRx X9:cloneJ'
rm(de_genes)
de_genes_535cx <- read.csv(paste0(output_dir,'signif_genes.csv'), check.names = F, stringsAsFactors = F)
de_genes <- de_genes_535cx
dim(de_genes)
summary(de_genes$logFC)
xl <- c(-3.5,max(de_genes$logFC))
topGenes <- c(as.character(markers_ls_upreg$gene_symbol),as.character(markers_ls_downreg$gene_symbol))


vol_535_cx <- plot_DE_genes_ggplot(de_genes, topGenes, capstr='', FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                                    plttitle, save_dir=output_dir, legendVisible=F,
                                    iscaption=TRUE, legend_verbose='none', save_plot=TRUE, xl=xl)
# vol_535_cx
# ggdraw(vol_535_cx$p_vol)

vol_535_cx_legend <- plot_DE_genes_ggplot(de_genes, topGenes, capstr='', FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                                          plttitle, save_dir=output_dir, legendVisible=T,
                                          iscaption=TRUE, legend_verbose='top', save_plot=TRUE, yl=c())





plotls <- list(vol_609$p_vol, vol_1035$p_vol)
plotls1 <- list(vol_535_cis$p_vol, vol_535_cx$p_vol)

volplotls <- cowplot::plot_grid(vol_609$p_vol, vol_535_cis$p_vol, vol_1035$p_vol, ncol=1, align = 'v')
prop_pls <- cowplot::plot_grid(vol_609$p_prop,vol_535_cis$p_prop,vol_1035$p_prop, ncol=1,
                               align = 'v')
prop_pls <- cowplot::plot_grid(vol_609$p_prop_cis,vol_609$p_prop_trans,
                               vol_535_cis$p_prop_cis,vol_535_cis$p_prop_trans,
                               vol_1035$p_prop_cis, vol_1035$p_prop_trans, ncol=1,
                               align = 'v')
plg <- cowplot::get_legend(vol_535_cis_legend)
p <- cowplot::plot_grid(volplotls, prop_pls, ncol=2,rel_widths = c(4,2.5)) #align = 'h', ,rel_widths = c(4,2)

ptotal <- cowplot::plot_grid(p, plg, ncol=1,rel_heights = c(10,1.5))
save_dir <- paste0(output_dir,'volcanos/')
dir.create(save_dir)
png(paste0(save_dir,"SUPP_Fig22_volcanos.png"), height = 2*850, width=2*650, res = 2*72)
print(ptotal)
dev.off()

ggsave(paste0(save_dir, "SUPP_Fig22_volcanos.pdf"),
       plot = ptotal,
       height = 12,
       width = 9.5,
       units = 'in',
       useDingbats=F)

ggsave(paste0(save_dir,"SUPP_Fig22_volcanos.png"),
       plot = ptotal,
       height = 9,
       width = 7.5,
       # useDingbats=F,
       type = "cairo-png",
       dpi=200
)


p <- cowplot::plot_grid(plotlist = plotls, ncol=2,  labels =  c('a','b'), align = 'hv') #align = 'h',
p1 <- cowplot::plot_grid(plotlist = plotls1, ncol=2,  labels =  c('c','d'), align = 'hv') #align = 'h',

plotls_lg <- list(vol_609$p_prop, NULL,NULL,vol_1035$p_prop)
p_lg <- cowplot::plot_grid(plotlist = plotls_lg, ncol=4, rel_widths = c(1,1.5,1.5,1))


plotls_lg1 <- list(vol_535_cis$p_prop, cowplot::ggdraw() + cowplot::draw_plot(lg),vol_535_cx$p_prop)
p_lg1 <- cowplot::plot_grid(plotlist = plotls_lg1, ncol=3, rel_widths =  c(1,3,1))


p_total <- cowplot::plot_grid(p, p_lg, p1, p_lg1, ncol=1, rel_heights = c(1, 0.2, 1, 0.2))

save_fig_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/trackplots_v3/"
# png(paste0(save_fig_dir, "Fig3_volcanoplots.png"), height =2*1050, width= 2*750,res = 2*72)
# print(p_total)
# dev.off()

# I found out that Farhia's thesis text width is 6.5 inches and text height is 9 inches. 
# We need that information to generate the plots. For example if I want to put 2 cloud 
# plots next to each other, I have to save  them as pdf about 3.2 inches wide each. 
# I managed to do this for the cloud plots, see below, and then gene label size is 6, 
# label axes and text are 8/10 etc. 
# I don't manage to fit as many labels as before, but it will have to do.
ggsave(paste0(save_fig_dir, "Fig5.9_volcanoplots_24_Feb.pdf"),
       plot = p_total,
       height = 10,
       width = 8.5,
       units = 'in',
       useDingbats=F)

ggsave(paste0(save_fig_dir,"Fig5.9_volcanoplots_24_Feb.png"),
       plot = p_total,
       height = 11,
       width = 9,
       # useDingbats=F,
       type = "cairo-png",
       dpi=200
)



plotls <- list(vol_609$p_vol, vol_1035$p_vol, vol_535_cis$p_vol)


p <- cowplot::plot_grid(plotlist = plotls, ncol=3,  labels =  c('a','b','c'), align = 'hv') #align = 'h',
# lg <- cowplot::get_legend(vol_535_cis_legend)

plotls_lg <- list(vol_609$p_prop, NULL,NULL,vol_1035$p_prop,NULL,NULL,vol_535_cis$p_prop)
p_lg <- cowplot::plot_grid(plotlist = plotls_lg, ncol=7, rel_widths = c(1,1.5,1.5,1,1.5,1.5,1))
p_total <- cowplot::plot_grid(p, p_lg, ncol=1, rel_heights = c(1, 0.2))
ggsave(paste0(save_dir,"volcanoplots.png"),
       plot = p_total,
       height = 7,
       width = 15,
       # useDingbats=F,
       type = "cairo-png",
       dpi=200
)




trans_5x <- de_genes %>%
  dplyr::filter(is.na(classified_gene_dlp))


# 0.58  0.35 vs 0.51  0.3
round(mean(abs(cis_5x$logFC)),2) # 0.58
round(sd(abs(cis_5x$logFC)),2)   # 0.35

round(mean(abs(trans_5x$logFC)),2) # 0.51
round(sd(abs(trans_5x$logFC)),2) # 0.3

stat_5x <- compute_stat(de_genes, 'SA535_CX5461')
stat_5x_nb <- compute_stat_nb(de_genes, 'SA535_CX5461')
stat_5x_top <- compute_stat_top50(de_genes, 'SA535_CX5461')
cis_5x_top <- compute_stat_top50(de_genes, 'in-cis')
top5pct_535cx <- compute_stat_top5percent(de_genes_535cx, 'SA535_CX5461')


incis_common1 <- intersect(top5pct_609$incis_genes,top5pct_1035$incis_genes)
incis_common2 <- intersect(top5pct_535cx$incis_genes, top5pct_535cis$incis_genes)
incis_common <- intersect(incis_common1, incis_common2)
incis_common <- intersect(top5pct_609$incis_genes, incis_common2)
intrans_common1 <- intersect(top5pct_609$intrans_genes,top5pct_1035$intrans_genes)
intrans_common2 <- intersect(top5pct_535cx$intrans_genes, top5pct_535cis$intrans_genes)
intrans_common <- intersect(intrans_common1,intrans_common2)

intrans_common <- intersect(top5pct_609$intrans_genes,intrans_common2)
# t <- t.test(abs(trans_5x$logFC), abs(cis_5x$logFC), alternative = 'greater')
# t$p.value < 0.05
stat <- as.data.frame(dplyr::bind_rows(stat_609, stat_1035, stat_5c, stat_5x))
View(stat)
write.csv(stat, paste0(base_dir, 'rnaseq_v6/trackplots_v3/stat_volcanos.csv'), quote=F, row.names = F)

stat_top50 <- as.data.frame(dplyr::bind_rows(stat_609_top, stat_1035_top, 
                                             stat_5c_top, stat_5x_top))
rownames(stat_top50) <- c('609_trans','609_cis','1035_trans','1035_cis',
                          '535c_trans','535c_cis','535x_trans','535x_cis')
View(stat_top50)
write.csv(stat, paste0(base_dir, 'rnaseq_v6/trackplots_v3/stat_volcanos_top50.csv'), quote=F, row.names = F)



compute_stat <- function(de_genes, datatag){
  cis <- de_genes %>%
    dplyr::filter(!is.na(classified_gene_dlp))
  
  trans <- de_genes %>%
    dplyr::filter(is.na(classified_gene_dlp))
  
  t <- t.test(abs(cis$logFC), abs(trans$logFC), alternative = 'greater')

  is_signif <- t$p.value < 0.05
  return(list(datatag=datatag, 
              cis_avg=round(mean(abs(cis$logFC)),2), cis_sd=round(sd(abs(cis$logFC)),2),
              trans_avg=round(mean(abs(trans$logFC)),2), trans_sd=round(sd(abs(trans$logFC)),2),
              is_signif=is_signif))
  
}

compute_stat_nb <- function(de_genes, datatag){
  cis <- de_genes %>%
    dplyr::filter(!is.na(classified_gene_dlp))
  
  trans <- de_genes %>%
    dplyr::filter(is.na(classified_gene_dlp))
  
  
  return(list(datatag=datatag, 
              nb_cis=nrow(cis),
              pct_cis=round(nrow(cis)/nrow(de_genes) * 100,2),
              nb_trans=nrow(trans),
              pct_trans=round(nrow(trans)/nrow(de_genes) * 100,2)))
  
}

compute_stat_top50 <- function(de_genes, gene_type){  #'in-trans'
  nbtopup=25
  nbtopdown=25
  minLogFC <- 0.25
  if(gene_type=='in-trans'){
    de_genes <- de_genes %>%
      dplyr::filter(is.na(classified_gene_dlp))
  }else{
    de_genes <- de_genes %>%
      dplyr::filter(!is.na(classified_gene_dlp))
  }
  
  markers_ls_upreg <- de_genes[de_genes$logFC>minLogFC,]
  markers_ls_upreg <- markers_ls_upreg[order(markers_ls_upreg$logFC,decreasing = T),] 
  markers_ls_upreg <- markers_ls_upreg[1:nbtopup,]
  markers_ls_downreg <- de_genes[de_genes$logFC<(-minLogFC),]
  markers_ls_downreg <- markers_ls_downreg[order(markers_ls_downreg$logFC,decreasing = F),] 
  markers_ls_downreg <- markers_ls_downreg[1:nbtopdown,]
  markers_ls <- dplyr::bind_rows(markers_ls_upreg, markers_ls_downreg)
  
  
  return(list(gene_type=gene_type, avg=round(mean(abs(markers_ls$logFC)),2), sd=round(sd(abs(markers_ls$logFC)),2)))
  
}

compute_stat_top5percent <- function(de_genes, datatag){
  
  markers_ls <- de_genes[order(abs(de_genes$logFC),decreasing = T),] 
  pct_top <- 0.05
  nbtop <- round(pct_top*nrow(markers_ls))
  print(nbtop)
  markers_ls <- markers_ls[1:nbtop,]
  cis <- markers_ls %>%
    dplyr::filter(!is.na(classified_gene_dlp))
  trans <- markers_ls %>%
    dplyr::filter(is.na(classified_gene_dlp))
  return(list(datatag=datatag, 
              pct_cis=round(nrow(cis)/nrow(markers_ls)*100,2),
              incis_genes=as.character(cis$gene_symbol),
              intrans_genes=as.character(trans$gene_symbol),
              pct_trans=round(nrow(trans)/nrow(markers_ls)*100,2)))
}

compute_stat_top50 <- function(de_genes, datatag){  #'in-trans'
    de_genes_t <- de_genes %>%
      dplyr::filter(is.na(classified_gene_dlp))
    mk_t <- get_topmarkers(de_genes_t)
    trans <- list(datatag=datatag, gene_type='in-trans', 
                  avg=round(mean(abs(mk_t$logFC)),2), sd=round(sd(abs(mk_t$logFC)),2))
    de_genes_c <- de_genes %>%
      dplyr::filter(!is.na(classified_gene_dlp))
    mk_c <- get_topmarkers(de_genes_c)
    cis <- list(datatag=datatag, gene_type='in-cis', 
                  avg=round(mean(abs(mk_c$logFC)),2), sd=round(sd(abs(mk_c$logFC)),2))
    stat_top50 <- as.data.frame(dplyr::bind_rows(trans, cis))
    
    t <- t.test(mk_t$logFC, mk_c$logFC, alternative = 'greater')
    stat_top50$is_significant <- t$p.value < 0.05
    return(stat_top50)
}
get_topmarkers <- function(de_genes){
  nbtopup=25
  nbtopdown=25
  minLogFC <- 0.25
  markers_ls_upreg <- de_genes[de_genes$logFC>minLogFC,]
  markers_ls_upreg <- markers_ls_upreg[order(markers_ls_upreg$logFC,decreasing = T),] 
  markers_ls_upreg <- markers_ls_upreg[1:nbtopup,]
  markers_ls_downreg <- de_genes[de_genes$logFC<(-minLogFC),]
  markers_ls_downreg <- markers_ls_downreg[order(markers_ls_downreg$logFC,decreasing = F),] 
  markers_ls_downreg <- markers_ls_downreg[1:nbtopdown,]
  markers_ls <- dplyr::bind_rows(markers_ls_upreg, markers_ls_downreg)
  return(markers_ls)
}



markers_ls_upreg <- de_genes[de_genes$logFC>minLogFC,]
markers_ls_upreg <- markers_ls_upreg[order(markers_ls_upreg$logFC,decreasing = T),] 
markers_ls_upreg <- markers_ls_upreg[1:nbtopup,]
markers_ls_downreg <- de_genes[de_genes$logFC<(-minLogFC),]
markers_ls_downreg <- markers_ls_downreg[order(markers_ls_downreg$logFC,decreasing = F),] 
markers_ls_downreg <- markers_ls_downreg[1:nbtopdown,]
topGenes <- c(as.character(markers_ls_upreg$gene_symbol),as.character(markers_ls_downreg$gene_symbol))


# vol_535_cx <- plot_DE_genes_edgeR(de_genes, topGenes, capstr='', FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
#                                   plttitle, save_dir=output_dir, legendVisible=T,
#                                   iscaption=TRUE, legend_verbose='right', save_plot=TRUE)
# 
# vol_535_cx

vol_535_cx_legend <- plot_DE_genes_edgeR(de_genes, topGenes, capstr='', 
                                         FDRcutoff=0.01, logFCcutoff=0.25,
                                         plttitle, output_dir, legendVisible=T,
                                         iscaption=TRUE)





lg <- cowplot::get_legend(vol_535_cx_legend)

pc <- ggdraw() + draw_plot(lg)
bottom_row <- plot_grid(pc, NULL, NULL, ncol = 3)




save_dir <- paste0(input_dir,'rnaseq_v6/normalization_evaluation/')
pv609 <- readRDS(paste0(save_dir,'SA609/edgeR_DE_analysis_SA609:_UTTTT-A_vs_UUUUU-H.rds'))
pv1035 <- readRDS(paste0(save_dir,'SA1035/edgeR_DE_analysis_SA609:_UTTTT-A_vs_UUUUU-H.rds'))
pv535_cis <- readRDS(paste0(save_dir,'SA535_cisplatin/edgeR_DE_analysis_SA535_cisplatin:_UUTTTT-T_vs_UUUUUU-J.rds'))
pv535_cx <- readRDS(paste0(save_dir,'SA535_CX5461/edgeR_DE_analysis_SA535_CX5461:_UXXXX-U_vs_UUUUU-J.rds'))
plotls <- list(pv609, pv1035, pv535_cis, pv535_cx)

lbs <- c('a','b','c','d')

p_total <- cowplot::plot_grid(plotlist = plotls, ncol=1,  align = 'v', labels = lbs)


png(paste0(save_dir, "eval_edgeR_using_scTransform.png"), height =2*1500, width= 2*1100,res = 2*72)
print(p_total)
dev.off()






