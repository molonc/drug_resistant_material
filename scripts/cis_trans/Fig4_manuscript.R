# 
# # main script: run_incis_intrans_analysis_v2.R
# 
# viz_Fig4 <- function(plot_ls, prop_plt_ls, pathway_plt_ls, save_dir){
#   
#   # dotplot here 
#   # save_dir <- input_dir
#   stat <- data.table::fread(paste0(save_dir,'summary_incis_intrans_genes_plots.csv')) %>% as.data.frame()
#   dim(stat)
#   # stat1 <- stat1 %>% left_join(gt, by='gene_type')
#   # stat1$pct_genes
#   # stat_backup <- stat1
#   # pd <- 'SA535'
#   unique(stat1$comp_type)
#   meta_ptx <- data.frame(datatag=c("SA501","SA530","SA604","SA609","SA535","SA1035"), 
#                          pt=paste0("Pt",rep(1:6,1)))
#   stat <- stat %>% left_join(meta_ptx, by='datatag')
#   dim(stat)
#   
#   # stat$plt_desc <- paste0(stat$pt,':\n',stat$plt_desc)
#   stat$plt_desc <- stat$labels_detailed
#   stat$desc <- stat$file_header
#   pdxs = c('SA609','SA535','SA1035')
#   pdxs_untreated = c('SA501','SA530','SA604','SA609','SA535','SA1035')
#   my_font <- "Helvetica"
#   
#   
#   # First get untreated comparison
#   df_untreated <- stat %>% 
#     dplyr::filter(datatag %in% pdxs_untreated & comp_type=='untreated')
#   # dim(df_untreated)
#   # unique(df_untreated$desc)
#   
#   # df_untreated$ord <- rep(1:dim(df_untreated)[1],1)
#       
#   plot_ls <- list()
#   
#   # df_untreated$plt_desc <- factor(df_untreated$plt_desc, levels = unique(df_untreated$plt_desc))
#   
#   # Untreated series ## TO DO: add patient names here
#   plot_ls[['untreated']] <- dotplot_prevalence(df_untreated)
#   df_untreated$desc <- df_untreated$file_header
#   res_untreated_plt <- plot_cis_pos_neg(df_untreated)
#   plot_ls[['untreated_proportion']] <- res_untreated_plt$p
#   
#   # Time series
#   for(pd in pdxs){
#     df <- stat %>% 
#       dplyr::filter(datatag==pd & comp_type=='treated_vs_untreated')
#     dim(df)
#         
#     plot_ls[[pd]] <- dotplot_prevalence(df)
#     res_plt <- plot_cis_pos_neg(df)
#     plot_ls[[paste0(pd,'_proportion')]] <- res_plt$p
#     
#   }
#   
#   p <- ggplot(df, aes(x = pct_genes, y = plt_desc, color = Gene_Type)) + 
#     geom_point(size=5, alpha=1) + #, size=2*log10(pct_genes)
#     scale_color_manual(name = NULL, values = keyvals_colour) + 
#     theme(legend.position ='bottom',
#           legend.text = element_text(size = 9, hjust = 0.5, family=my_font),
#           legend.title = element_text(size = 10, hjust = 0.5, family=my_font, angle = 90))
#   p <- p + guides(color = guide_legend(override.aes = list(size=8, ncol = 2)))#shape = 0, , ncol = 2
#   lg <- cowplot::get_legend(p)
#   plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
#   # plg
#   
#   
#   
#   plot_ls[['lg']] <- plg
#   plot_ls[['lg_cis']] <- res_untreated_plt$plg
#   
#   # lgs <- cowplot::plot_grid(prop_plt$plg, pathway_plt$plg, ncol=1)
#   # lgs2 <- cowplot::plot_grid(plot_ls$lg, plot_ls$lg_cis, ncol=1)
#   lgs <- cowplot::plot_grid(prop_plt_ls$plg, pathway_plt_ls$plg, ncol=1)
#   lgs2 <- cowplot::plot_grid(plot_ls$lg, plot_ls$lg_cis, ncol=1)
#   lgs_total <- cowplot::plot_grid(lgs2, lgs, nrow = 2)
#   untreated_row <- cowplot::plot_grid(plot_ls$untreated, plot_ls$untreated_proportion, 
#                                       lgs_total,
#                                       rel_widths = c(2.8,1,2.5), ncol=3)
#   sa609_row <- cowplot::plot_grid(plot_ls$SA609, plot_ls$SA609_proportion, 
#                                   prop_plt_ls$SA609, pathway_plt_ls$SA609,
#                                   rel_widths = c(2.8,1,2,0.5), ncol=4)
#   sa535_row <- cowplot::plot_grid(plot_ls$SA535, plot_ls$SA535_proportion, 
#                                   prop_plt_ls$SA535, pathway_plt_ls$SA535,
#                                   rel_widths = c(2.8,1,2,0.5), ncol=4)
#   sa1035_row <- cowplot::plot_grid(plot_ls$SA1035, plot_ls$SA1035_proportion, 
#                                    prop_plt_ls$SA1035, pathway_plt_ls$SA1035,
#                                    rel_widths = c(2.8,1,2,0.5), ncol=4)
#   # lgt_row <- cowplot::plot_grid(plot_ls$lg, plot_ls$lg_cis, rel_widths = c(4,1,2,1), ncol=4)
#   # p_total <- cowplot::plot_grid(untreated_row,sa609_row,sa535_row,sa1035_row,lgt_row,
#   #                               ncol = 1,align='vh', rel_heights = c(6,4,3,2,1))
#   p_total <- cowplot::plot_grid(NULL, untreated_row,NULL,sa609_row,NULL,sa535_row,NULL,sa1035_row,
#                                 ncol = 1,align='v', rel_heights = c(0.35,6,0.35,7,0.35,4,0.35,5),
#                                 labels = c('a -Rx PDX tumors Pt1-Pt6','',
#                                            'b -Rx/+Rx time series Pt4','',
#                                            'c -Rx/+Rx time series Pt5','',
#                                            'd -Rx/+Rx time series Pt6',''))
#   save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/Fig4_prevalence_cistrans/'
#   # dir.create(save_fig_dir)
#   saveRDS(untreated_row, paste0(save_fig_dir,'untreated_row.rds'))
#   saveRDS(sa609_row, paste0(save_fig_dir,'sa609_row.rds'))
#   saveRDS(sa535_row, paste0(save_fig_dir,'sa535_row.rds'))
#   saveRDS(sa1035_row, paste0(save_fig_dir,'sa1035_row.rds'))
#   # saveRDS(lgt_row, paste0(save_fig_dir,'lgt_row.rds'))
#   saveRDS(p_total, paste0(save_fig_dir,'p_total.rds'))
#   
#   # p_total <- cowplot::plot_grid(plotlist=plot_ls, ncol = 1, 
#   #                               align='vh', rel_heights = c(6,4,3,2,1))
#   # labels = c('UnRx vs. UnRx','SA609: Res Rx vs. Sen UnRx',
#   #            'SA535: Res Rx vs. Sen UnRx','SA1035: Res Rx vs. Sen UnRx')
#   # labels = c('UnRx vs. UnRx','Rx vs. UnRx','Rx vs. UnRx','Rx vs. UnRx'),
#   # label_size = 10,
#   # label_fontfamily = my_font,
#   # label_x = 0, label_y = 0
#   # )#hjust = 0, vjust = 0
#   # p_total <- p_total + labs(x=NULL, y="(%) Gene Type ", title='Differentially Expressed Genes - Copy Number Driven Proportion')
#   # png(paste0(save_fig_dir,"Fig4_incis_intrans_genes_prevalence_summary.png"), 
#   #     height = 2*800, width=2*750, res = 2*72)
#   # print(p_total)
#   # dev.off()
#   
#   ggsave(paste0(save_fig_dir,"Fig4_incis_intrans_genes_prevalence_summary.pdf"),
#          plot = p_total,
#          height = 14,
#          width = 9.8,
#          useDingbats=F,
#          dpi=150)
#   ggsave(paste0(save_fig_dir,"Fig4_incis_intrans_genes_prevalence_summary.png"),
#          plot = p_total,
#          height = 14,
#          width = 9.8,
#          # useDingbats=F,
#          dpi=150)
#   # saveRDS(p_total, paste0(save_dir,"incis_intrans_genes_prevalence_dotplot.rds"))
# }