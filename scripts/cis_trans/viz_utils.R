library(viridis)
library(dplyr)


my_font <- "Helvetica"

ct_theme <-   
  theme_bw(base_size=11) + 
  theme(#legend.position = "bottom",
    legend.position = "none",
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor.x = element_blank(),
    # axis.ticks.y = element_blank(),
    # axis.text  = element_blank(),
    text = element_text(size = 11, hjust = 0.5, family=my_font),
    # axis.text.x = element_text(size=11, hjust = 0.5, family=my_font, angle = 90),  #
    axis.text.x = element_blank(),
    # axis.text.x = element_text(size=11, hjust = 0.5, family=my_font, angle = 90, color='black'), ## for testing
    # axis.text.y = element_text(size=11, hjust = 0.5, family=my_font), ## for testing
    axis.text.y = element_blank(), #element_text(size=11, hjust=0.5, family=my_font),
    plot.title = element_text(size=11, hjust=0.5, family=my_font), #, face="bold"
    axis.title.y = element_blank(),
    # axis.title.x = element_text(size=11, hjust = 0.5, family=my_font)
    # axis.title.x = element_blank()
    # strip.text.x = element_text(color="black",size=11, family=my_font),
    # strip.text.y = element_text(color="black",size=11, family=my_font),
    strip.text.x = element_blank(), 
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    # strip.placement = "outside",
    panel.background = element_rect(fill = "white", colour = "grey92"),
    # panel.grid = element_line(colour = "grey92"), 
    panel.grid.major = element_line(colour = "grey92"),
    # panel.grid.minor = element_line(colour = "grey92", size = rel(0.5)),
    # strip.background = element_rect(fill = "grey85", colour = "grey20") ,
    panel.border = element_rect(fill = NA, 
                                colour = "grey92")
  ) #+

viz_Fig4_v2 <- function(plot_ls, prop_plt_ls, pathway_plt_ls, save_dir){
  rel_ht <- c(5,2,5)
  rx_lg <- cowplot::plot_grid(plot_ls$treated_SA609_c1lb, plot_ls$treated_SA535_c1lb,
                              plot_ls$treated_SA1035_c1lb, rel_heights = rel_ht, ncol=1)
  unrx_lg <- cowplot::plot_grid(plot_ls$treated_SA609_c2lb, plot_ls$treated_SA535_c2lb,
                                plot_ls$treated_SA1035_c2lb, rel_heights = rel_ht, ncol=1)
  
  passage_lg <- cowplot::plot_grid(plot_ls$treated_SA609_passage, plot_ls$treated_SA535_passage,
                                   plot_ls$treated_SA1035_passage, rel_heights = rel_ht, ncol=1)
  
  # prop_plt_ls, pathway_plt_ls
  plt_ct <- cowplot::plot_grid(passage_lg, rx_lg, unrx_lg,
                               plot_ls$treated,plot_ls$treated_proportion,plot_ls$treated_cistrans,
                               prop_plt_ls$prevalence_plt, pathway_plt_ls$pathway_sig_plt,
                               rel_widths = c(0.15,0.15,0.15,
                                              1,1,1,
                                              1,0.3), nrow=1)
  ggsave(paste0(save_dir,"Fig23_cistrans.png"),
         plot = plt_ct,
         height = 5,
         width = 8,
         # useDingbats=F,
         dpi=300)
  return(plt_ct)
  
  
}

viz_SUPP_Fig6_DrugHolidayRxH <- function(plot_ls, plot_untreated_ls, prop_plt_ls, pathway_plt_ls, save_dir){
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  untreated_passage_lg <- cowplot::plot_grid(plot_untreated_ls$UnRx_SA501_passage, plot_untreated_ls$UnRx_SA530_passage,
                                           plot_untreated_ls$UnRx_SA604_passage, plot_untreated_ls$UnRx_SA609_passage,
                                           plot_untreated_ls$UnRx_SA535_passage, plot_untreated_ls$UnRx_SA1035_passage,
                                           ncol=1)
  
  untreated_unrx1_lg <- cowplot::plot_grid(plot_untreated_ls$UnRx_SA501_c1lb, plot_untreated_ls$UnRx_SA530_c1lb,
                                             plot_untreated_ls$UnRx_SA604_c1lb, plot_untreated_ls$UnRx_SA609_c1lb,
                                             plot_untreated_ls$UnRx_SA535_c1lb, plot_untreated_ls$UnRx_SA1035_c1lb,
                                             ncol=1)
  untreated_unrx2_lg <- cowplot::plot_grid(plot_untreated_ls$UnRx_SA501_c2lb, plot_untreated_ls$UnRx_SA530_c2lb,
                                           plot_untreated_ls$UnRx_SA604_c2lb, plot_untreated_ls$UnRx_SA609_c2lb,
                                           plot_untreated_ls$UnRx_SA535_c2lb, plot_untreated_ls$UnRx_SA1035_c2lb,
                                           ncol=1)
  plg_total <- cowplot::plot_grid(plot_untreated_ls$lg_cis, plot_untreated_ls$lg_cistrans, nrow=2)
  untreated_plt_ct <- cowplot::plot_grid(untreated_passage_lg, untreated_unrx1_lg, untreated_unrx2_lg,
                               plot_untreated_ls$UnRx,plot_untreated_ls$UnRx_proportion,plot_untreated_ls$UnRx_cistrans,
                               plg_total,
                               rel_widths = c(0.15,0.15,0.15,
                                              1,1,1,1), nrow=1)
  
  
  
  
  rel_ht <- c(5,2,5)
  rx_lg <- cowplot::plot_grid(plot_ls$treated_SA609_c1lb, plot_ls$treated_SA535_c1lb,
                              plot_ls$treated_SA1035_c1lb, rel_heights = rel_ht, ncol=1)
  unrx_lg <- cowplot::plot_grid(plot_ls$treated_SA609_c2lb, plot_ls$treated_SA535_c2lb,
                                plot_ls$treated_SA1035_c2lb, rel_heights = rel_ht, ncol=1)
  
  passage_lg <- cowplot::plot_grid(plot_ls$treated_SA609_passage, plot_ls$treated_SA535_passage,
                                   plot_ls$treated_SA1035_passage, rel_heights = rel_ht, ncol=1)
  
  # prop_plt_ls, pathway_plt_ls
  
  plt_ct <- cowplot::plot_grid(passage_lg, rx_lg, unrx_lg,
                               plot_ls$treated,plot_ls$treated_proportion,plot_ls$treated_cistrans,
                               prop_plt_ls$prevalence_plt, pathway_plt_ls$pathway_sig_plt,
                               rel_widths = c(0.15,0.15,0.15,
                                              1,1,1,1), nrow=1) #,0.3
  
  plt_SUPPFig61 <- cowplot::plot_grid(NULL, untreated_plt_ct, NULL, plt_ct, NULL, ncol=1, rel_heights=c(1,1,0.2,2,0.2))
  plt_SUPPFig6 <- cowplot::plot_grid(NULL, plt_SUPPFig61, NULL, nrow=1, rel_widths = c(0.5,5, 0.1))
  
  ggsave(paste0(save_dir,"SUPP_Fig6_cistrans_RxH_UnRx.svg"),
         plot = plt_SUPPFig6,
         height = 11,
         width = 9,
         # useDingbats=F,
         dpi=250)
  
  return(plt_ct)
  
  
}
viz_Fig4 <- function(plot_ls, prop_plt_ls, pathway_plt_ls, save_dir){
  save_fig_dir <- paste0(save_dir, 'figs/')
  if(!dir.exists(save_fig_dir)){
    dir.create(save_fig_dir)
  }
  
  lgs <- cowplot::plot_grid(prop_plt_ls$plg, pathway_plt_ls$plg, ncol=1)
  lgs2 <- cowplot::plot_grid(plot_ls$lg, plot_ls$lg_cis, ncol=1)
  
  lgs_total <- cowplot::plot_grid(lgs2, lgs, nrow = 2)
  
  ## version 1
  # untreated_row <- cowplot::plot_grid(plot_ls$untreated, plot_ls$untreated_proportion, 
  #                                     lgs_total,
  #                                     rel_widths = c(2.8,1.2,2.5), ncol=3)
  # sa609_row <- cowplot::plot_grid(plot_ls$SA609, plot_ls$SA609_proportion, 
  #                                 prop_plt_ls$SA609, pathway_plt_ls$SA609,
  #                                 rel_widths = c(2.8,1.2,2,0.5), ncol=4)
  # sa535_row <- cowplot::plot_grid(plot_ls$SA535, plot_ls$SA535_proportion, 
  #                                 prop_plt_ls$SA535, pathway_plt_ls$SA535,
  #                                 rel_widths = c(2.8,1.2,2,0.5), ncol=4)
  # sa1035_row <- cowplot::plot_grid(plot_ls$SA1035, plot_ls$SA1035_proportion, 
  #                                  prop_plt_ls$SA1035, pathway_plt_ls$SA1035,
  #                                  rel_widths = c(2.8,1.2,2,0.5), ncol=4)
  
  
  untreated_row <- cowplot::plot_grid(plot_ls$untreated_passage,plot_ls$untreated_c1lb,plot_ls$untreated_c2lb,
                                      plot_ls$untreated, plot_ls$untreated_proportion, 
                                      NULL,
                                      rel_widths = c(0.3,0.3,0.3,2,1.2,3.1), ncol=6)
  sa609_row <- cowplot::plot_grid(plot_ls$SA609_passage,plot_ls$SA609_c1lb,plot_ls$SA609_c2lb,
                                  plot_ls$SA609, plot_ls$SA609_proportion, plot_ls$SA609_cistrans,
                                  prop_plt_ls$SA609, pathway_plt_ls$SA609,
                                  rel_widths = c(0.3,0.3,0.3,2,1.2,1.2,1.4,0.5), ncol=8)
  sa535_row <- cowplot::plot_grid(plot_ls$SA535_passage,plot_ls$SA535_c1lb,plot_ls$SA535_c2lb,
                                  plot_ls$SA535, plot_ls$SA535_proportion, plot_ls$SA535_cistrans,
                                  prop_plt_ls$SA535, pathway_plt_ls$SA535,
                                  rel_widths = c(0.3,0.3,0.3,2,1.2,1.2,1.4,0.5), ncol=8)
  sa1035_row <- cowplot::plot_grid(plot_ls$SA1035_passage,plot_ls$SA1035_c1lb,plot_ls$SA1035_c2lb,
                                   plot_ls$SA1035, plot_ls$SA1035_proportion, plot_ls$SA1035_cistrans,
                                   prop_plt_ls$SA1035, pathway_plt_ls$SA1035,
                                   rel_widths = c(0.3,0.3,0.3,2,1.2,1.2,1.4,0.5), ncol=8)
  
  
  # lgt_row <- cowplot::plot_grid(plot_ls$lg, plot_ls$lg_cis, rel_widths = c(4,1,2,1), ncol=4)
  # p_total <- cowplot::plot_grid(untreated_row,sa609_row,sa535_row,sa1035_row,lgt_row,
  #                               ncol = 1,align='vh', rel_heights = c(6,4,3,2,1))
  
  # theme_set(cowplot::theme_cowplot(font_size=11, font_family = "Helvetica"))
  p_total <- cowplot::plot_grid(NULL, untreated_row,NULL,
                                sa609_row,NULL,sa535_row,NULL,sa1035_row,
                                ncol = 1, align='v', rel_heights = c(0.5,6,0.5,
                                                                    7,0.5,4,0.5,5),
                                labels = c('a UnRx PDX tumors Pt1-Pt6','',
                                           'b UnRx/Rx time series Pt4','',
                                           'c UnRx/Rx time series Pt5','',
                                           'd UnRx/Rx time series Pt6',''))+
    theme(plot.background = element_rect(fill = "white", colour = "white")) 
    
  
  # png(paste0(save_fig_dir,'Fig4.png'), height = 2*1300, width=2*1600,res = 2*72)
  # print(p_total)
  # dev.off()
  
  ggsave(paste0(save_fig_dir,"Fig4_cistrans.png"),
         plot = p_total,
         height = 11,
         width = 12,
         # useDingbats=F,
         dpi=250)
  
  
  ggsave(paste0(save_fig_dir,"Fig4_cistrans.pdf"),
         plot = p_total,
         height = 13,
         width = 11,
         useDingbats=F,
         dpi=150)
  
  return(p_total)
  
}
viz_heatmap <- function(input_dir, output_dir){
  library(ComplexHeatmap)
  SA609_summary <- read.csv(paste0(input_dir,'SA609/SA609_result/SA609_summary.csv'), 
                            stringsAsFactors=F, check.names = F)
  # dim(SA609_summary)
  # View(SA609_summary)
  SA1035_summary <- read.csv(paste0(input_dir,'SA1035/SA1035_result/SA1035_summary.csv'), 
                             stringsAsFactors=F, check.names = F)
  
  SA535_cis_summary <- read.csv(paste0(input_dir,'SA535/SA535_cis_result/SA535_cis_summary.csv'), 
                                stringsAsFactors=F, check.names = F)
  SA535_cx_summary <- read.csv(paste0(input_dir,'SA535/SA535_cx_result/SA535_cx_summary.csv'), 
                               stringsAsFactors=F, check.names = F)
  
  
  stat_df <- list(SA609_summary, SA1035_summary, SA535_cis_summary, SA535_cx_summary)
  stat_df2 <- list()
  c <- 0
  for(df in stat_df){
    rownames(df) <- df$gene_type
    df <- df[,colnames(df) !='gene_type']
    exclude_cols <- grep('_nbgenes',colnames(df), value = T)
    cols_use <- colnames(df)[!colnames(df) %in% exclude_cols]
    df <- df[,cols_use]
    c <- c + 1
    stat_df2[[c]] <- t(df)
  }
  # View(SA535_cx_summary)
  stat <- do.call(rbind, stat_df2)
  # View(stat)
  # rownames(stat)
  coln <- colnames(stat)
  coln <- gsub('In_cis','InCis',coln)
  coln <- gsub('In_trans','InTrans',coln)
  colnames(stat) <- coln
  gt <- lapply(strsplit(colnames(stat), "_"), function(x) {
    return(x[1])
  })
  
  annotation_col = data.frame(
    Gene_Type = factor(as.character(gt))#, 
    #GeneType = 1:6
  )
  # dim(annotation_col)
  rownames(annotation_col) = colnames(stat)
  
  series <- lapply(strsplit(rownames(stat), "-"), function(x) {
    return(x[1])
  })
  annotation_row = data.frame(
    PDX = factor(as.character(series))
  )
  # dim(annotation_row)
  rownames(annotation_row) = rownames(stat)
  
  top_anno = ComplexHeatmap::HeatmapAnnotation(Gene_Type = factor(as.character(gt)),
                                               col = list(Gene_Type = c(InCis = "#FF3232", InTrans = "#40A0E0")))
  left_anno = ComplexHeatmap::rowAnnotation(PDX = factor(as.character(series)),
                                            col = list(PDX=c(SA609="#1A6A0B",SA1035 = "#8DD87F", SA535_cis = "#30F50A", SA535_cx = "#560C56"))  #SA609="#1A6A0B"
  )
  
  
  cell_func = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", stat[i, j]), x, y, gp = gpar(fontsize = 13))
  }
  
  # ha = HeatmapAnnotation(summary = anno_summary(height = unit(4, "cm")))
  p <- ComplexHeatmap::Heatmap(stat, na_col = "white",
                               show_column_names=T,
                               show_row_names = T,
                               cluster_rows=F,cluster_columns=F,
                               name = "% Genes", 
                               row_order = rownames(stat), ## TO DO: check rows order
                               row_split= annotation_row,
                               row_title_rot = 0,
                               row_gap = unit(2, "mm"),
                               column_split = annotation_col, 
                               column_title = "Genes Copy Number Driven",
                               column_gap = unit(2, "mm"),
                               column_names_gp = grid::gpar(fontsize = 10),
                               row_names_gp = grid::gpar(fontsize = 10),
                               show_heatmap_legend = T,
                               top_annotation=top_anno,
                               left_annotation = left_anno,
                               cell_fun = cell_func,
                               row_dend_reorder=F)
  
  
  png(paste0(output_dir,'summary_incis_intrans_genes.png'), height = 2*700, width=2*1000,res = 2*72)
  print(p)
  dev.off()
  
  # View(stat)
  write.csv(stat, paste0(output_dir,'summary_incis_intrans_genes.csv'), quote = F, row.names = T)
  
}

viz_prevalence <- function(input_dir, output_dir){
  stat <- read.csv(paste0(output_dir,'summary_incis_intrans_genes.csv'), check.names = F, stringsAsFactors=F, row.names = 1)
  # dim(stat)
  # View(stat)
  
  # stat$InCis_Increase <- stat$InCis_Increase_DownRegulated + stat$InCis_Increase_UpRegulated
  # stat$InCis_Decrease <- stat$InCis_Decrease_DownRegulated + stat$InCis_Decrease_UpRegulated
  # stat$InTrans <- stat$InTrans_DownRegulated + stat$InTrans_UpRegulated
  # in-cis up (proportion of CN change up)
  # in-cis down (proportion of CN change down)
  # in-trans (proportion of in-trans)
  # stat <- stat[,c('InCis_Increase','InCis_Decrease','InTrans')]
  
  stat1 <- stat
  # stat1$desc <- as.character(rownames(stat1))
  
  stat1 <- stat1 %>%
    pivot_longer(!desc, names_to = "gene_type", values_to = "percentage")
  
  series <- mclapply(strsplit(stat1$desc, "-"), function(x) {
    return(x[1])
  }, mc.cores = 2)
  
  desc <- mclapply(strsplit(stat1$desc, "-"), function(x) {
    return(x[2])
  }, mc.cores = 2)
  
  stat1$desc <- as.character(desc)
  stat1$series <- as.character(series)
  stat1$series <- factor(stat1$series, levels = unique(stat1$series))
  p <- ggplot(stat1, aes(fill=gene_type, y=percentage, x=desc)) + 
    geom_bar(position="fill", stat="identity")+
    facet_grid(. ~ series, scales="free", space='free')
  p <- p + labs(x='Resistant vs Sensitive cells comparisons', y="Percentage gene type (%)", title='Genes - Copy Number Driven Percentage')
  p <- p + theme(legend.title = element_text(size=10), 
                 legend.text = element_text(size=7),
                 plot.title = element_text(color="black", size=12, hjust = 0.5),
                 # legend.position = "none", 
                 # axis.line = element_blank(), 
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.border = element_blank(),
                 axis.text.x = element_text(size=10, angle = 90, color="black"),
                 axis.text.y = element_text(size=10, color="black")
                 #axis.title = element_blank(),
                 # axis.ticks = element_blank()
  )
  
  png(paste(output_dir,"incis_intrans_genes_prevalence_summary.png",sep="/"), height = 2*380, width=2*700, res = 2*72)
  print(p)
  dev.off()
}


plot_DE_genes_edgeR <- function(df, topGenes, capstr='', FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                                plttitle="A versus B", save_dir="",legendVisible=F,
                                iscaption=TRUE, legend_verbose='none', save_plot=TRUE){  
  #legend_verbose is none or 'right', 'left','bottom'
  # df <- de_genes
  # library(EnhancedVolcano)
  
  # colnames(df)[which(names(df) == "avg_logFC")] <- "log2FoldChange"
  # colnames(df)[which(names(df) == "p_val_adj")] <- "padj"
  # capstr <- paste0("FDR cutoff, ",FDRcutoff,"; logFC cutoff, ",logFCcutoff, "; nb genes signif, ",nbgenessig)
  # summary(as.factor(df$gene_type))
  Gene_Type=c('In_cis_Decrease_DownRegulated','In_cis_Decrease_UpRegulated',
              'In_cis_Increase_DownRegulated','In_cis_Increase_UpRegulated',
              'In_trans_DownRegulated','In_trans_UpRegulated'
  )
  gt <- data.frame(Gene_Type=Gene_Type,
                   gene_type=c('In-cis_DD','In-cis_DU',
                               'In-cis_ID','In-cis_IU',
                               'In-trans_D','In-trans_U'),
                   col_gt=c('#E61367','#D35B53',
                            '#EF816E','#E82503',
                            '#7BEADF','#093CB4'))
  
  df <- df %>% left_join(gt, by='Gene_Type')
  # unique(df$col_gt)
  keyvals_colour <- as.character(df$col_gt)
  names(keyvals_colour) <- df$gene_type
  # keyvals_colour <- factor(keyvals_colour, levels = unique(keyvals_colour))
  # unique(keyvals_colour)
  # names(keyvals_colour[1:3])
  df <- df[abs(df$logFC)>logFCcutoff & df$FDR<FDRcutoff  & df$PValue<pValuecutoff,]
  df$logFC <- sapply(df$logFC, function(x) replace(x, x > 3, 3))
  df$logFC <- sapply(df$logFC, function(x) replace(x, x < (-3), -3))
  
  st <- summary(df$gene_type)
  keyvals_shape <- ifelse(df$is_fitness_gene==T, 3, 16)
  names(keyvals_shape) <- ifelse(df$is_fitness_gene==T,'Fitness genes','Others')
  
  # df$gene_type <- factor(df$gene_type, levels = unique(df$gene_type))
  if(capstr=='' & iscaption){
    # capstr <- paste0(capstr,'With abs(logFC)>0.25, FDR<0.01, pValue<0.05 \n')
    capstr <- paste0(capstr,names(st[1]),':',as.numeric(st[1]), ', ')
    capstr <- paste0(capstr,names(st[2]),':',as.numeric(st[2]), ', ')
    capstr <- paste0(capstr,names(st[3]),':',as.numeric(st[3]), ', ')
    capstr <- paste0(capstr,names(st[4]),':',as.numeric(st[4]), ' \n')
    capstr <- paste0(capstr,names(st[5]),':',as.numeric(st[5]), ', ')
    capstr <- paste0(capstr,names(st[6]),':',as.numeric(st[6]), ' ')
    # for(i in rep(1:length(st),1)){
    #   # print(st[i])
    #   capstr <- paste0(capstr,names(st[i]),':',as.numeric(st[i]), ' ')
    # }
    
  }
  df$mlog10FDR <- -log10(df$FDR)
  
  df$mlog10FDR <- sapply(df$mlog10FDR, function(x) replace(x, is.infinite(x), 300))
  # legend_verbose <- 'top'
  p <- EnhancedVolcano::EnhancedVolcano(df,
                                        lab = df$gene_symbol,   #df$symbol or rownames(df)
                                        x = 'logFC',
                                        y = 'FDR',
                                        # selectLab = as.character(topGenes),
                                        selectLab = topGenes,
                                        xlim=c(-3,3),
                                        ylim = c(0, max(df$mlog10FDR, na.rm = TRUE)+10),
                                        xlab = bquote(~Log[2] ~ ' fold change'),
                                        ylab = bquote(~-Log[10]~italic(FDR)),
                                        pCutoff = FDRcutoff,
                                        FCcutoff = logFCcutoff,
                                        # pointSize = 3.2,
                                        pointSize = c(ifelse(df$is_fitness_gene==T, 4, 3)),
                                        labSize = 3,
                                        colAlpha = 0.8,
                                        gridlines.major = FALSE,
                                        gridlines.minor = FALSE,
                                        colCustom = keyvals_colour,
                                        shapeCustom = keyvals_shape,
                                        cutoffLineType = 'twodash',
                                        cutoffLineCol = 'grey0',
                                        cutoffLineWidth = 0.5,
                                        # hline = c(0.01),
                                        # hlineCol = c('grey0'),
                                        # hlineType = 'longdash',
                                        legend=NULL,
                                        # legend=c('NS','Log2 FC','Adjusted p-value',
                                        #          'Adjusted p-value & Log2 FC'),
                                        legendPosition = legend_verbose,
                                        legendLabSize = 8,
                                        legendIconSize = 3.0,
                                        legendVisible = legendVisible,
                                        # col=c('black', 'black', 'black', 'red3'),
                                        col=unique(keyvals_colour),
                                        title = NULL,
                                        subtitle = plttitle,
                                        caption = capstr,
                                        drawConnectors = TRUE,
                                        widthConnectors = 0.15,
                                        colConnectors = 'grey30')
  # p
  # p <- p + theme(legend.position="none")
  
  tag <- ifelse(legend_verbose=='none','','_with_legend')
  png(paste0(save_dir,"DE_",plttitle,tag,".png"), height = 2*550, width=2*650,res = 2*72)
  print(p)
  dev.off()
  if(save_plot){
    plttitle <- gsub(':','_',plttitle)
    plttitle <- gsub(' ','_',plttitle)
    saveRDS(p, file=paste0(save_dir,"DE_",plttitle,".rds"))
  }
  
  return(p)
}


plot_DE <- function(de, de_genes, pair_groups, output_dir,
                    minLogFC=0.25, pValueThrs=0.05,
                    nbtopup=30, nbtopdown=30){
  FDR_cutoff <- 0.01
  nbtop_extract <- 20
  de_genes <- de_genes[abs(de_genes$logFC)>minLogFC & de_genes$FDR<FDR_cutoff  & de_genes$PValue<pValueThrs,]
  # Plot DE genes
  # plttitle <- paste0(pair_groups[de,'datatag'],":  ",pair_groups[de,'clone1']," versus ",pair_groups[de,'clone2'])
  plttitle <- de
  
  markers_ls_upreg <- de_genes[de_genes$logFC>minLogFC,]
  markers_ls_upreg <- markers_ls_upreg[order(markers_ls_upreg$logFC,decreasing = T),] 
  dim(markers_ls_upreg)
  # pair_groups[de,'resistant_genes'] <- nrow(markers_ls_upreg)
  
  # Extract top incis, intrans up-regulated
  markers_ls_upreg_incis <- markers_ls_upreg[!is.na(markers_ls_upreg$classified_gene_dlp),]
  dim(markers_ls_upreg_incis)
  markers_ls_upreg_incis <- markers_ls_upreg_incis[order(markers_ls_upreg_incis$logFC,decreasing = T),] 
  if(nrow(markers_ls_upreg_incis) > nbtop_extract){
    markers_ls_upreg_incis <- markers_ls_upreg_incis[1:nbtop_extract,]
  }
  write.csv(markers_ls_upreg_incis, file=paste0(output_dir,'topgenes_upregulated_incis.csv'), quote = F, row.names = F)
  
  markers_ls_upreg_intrans <- markers_ls_upreg[is.na(markers_ls_upreg$classified_gene_dlp),]
  dim(markers_ls_upreg_intrans)
  markers_ls_upreg_intrans <- markers_ls_upreg_intrans[order(markers_ls_upreg_intrans$logFC,decreasing = T),] 
  if(nrow(markers_ls_upreg_intrans) > nbtop_extract){
    markers_ls_upreg_intrans <- markers_ls_upreg_intrans[1:nbtop_extract,]
  }
  write.csv(markers_ls_upreg_intrans, file=paste0(output_dir,'topgenes_upregulated_intrans.csv'), quote = F, row.names = F)
  
  
  # dim(markers_ls_upreg)
  if(nrow(markers_ls_upreg) < nbtopup){
    nbtopup <- nrow(markers_ls_upreg)
  }
  markers_ls_upreg <- markers_ls_upreg[1:nbtopup,]
  
  # Get top down-regulated genes
  markers_ls_downreg <- de_genes[de_genes$logFC<(-minLogFC),]
  markers_ls_downreg <- markers_ls_downreg[order(markers_ls_downreg$logFC,decreasing = F),] 
  
  # Extract top incis, intrans down-regulated
  markers_ls_downreg_incis <- markers_ls_downreg[!is.na(markers_ls_downreg$classified_gene_dlp),]
  dim(markers_ls_downreg_incis)
  markers_ls_downreg_incis <- markers_ls_downreg_incis[order(markers_ls_downreg_incis$logFC,decreasing = F),] 
  if(nrow(markers_ls_downreg_incis) > nbtop_extract){
    markers_ls_downreg_incis <- markers_ls_downreg_incis[1:nbtop_extract,]
  }
  write.csv(markers_ls_downreg_incis, file=paste0(output_dir,'topgenes_downregulated_incis.csv'), quote = F, row.names = F)
  
  markers_ls_downreg_intrans <- markers_ls_downreg[is.na(markers_ls_downreg$classified_gene_dlp),]
  dim(markers_ls_downreg_intrans)
  markers_ls_downreg_intrans <- markers_ls_downreg_intrans[order(markers_ls_downreg_intrans$logFC,decreasing = F),] 
  if(nrow(markers_ls_downreg_intrans) > nbtop_extract){
    markers_ls_downreg_intrans <- markers_ls_downreg_intrans[1:nbtop_extract,]
  }
  write.csv(markers_ls_downreg_intrans, file=paste0(output_dir,'topgenes_downregulated_intrans.csv'), quote = F, row.names = F)
  
  
  if(nrow(markers_ls_downreg) < nbtopdown){
    nbtopdown <- nrow(markers_ls_downreg)
  }
  markers_ls_downreg <- markers_ls_downreg[1:nbtopdown,]
  topGenes <- c(as.character(markers_ls_upreg$gene_symbol),as.character(markers_ls_downreg$gene_symbol))
  genes_df <- data.frame(desc=rep(de, length(topGenes)), topGenes=topGenes)
  write.csv(genes_df, file=paste0(output_dir,'topgenes.csv'), quote = F, row.names = F)
  
  # rownames(markers_ls_tmp) <- markers_ls_tmp$gene_symb
  
  plot_DE_genes_edgeR(de_genes, topGenes, capstr='', 
                      FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                      plttitle, output_dir, legendVisible=F,
                      iscaption=TRUE, legend_verbose='none', save_plot=TRUE)
  plot_DE_genes_edgeR(de_genes, topGenes, capstr='', 
                      FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                      plttitle, output_dir, legendVisible=T,
                      iscaption=TRUE, legend_verbose='top', save_plot=FALSE)
  
  
}

plot_dlp_prevalence <- function(output_dir=NULL){
  output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/SA535-v6/'
  base_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/')
  dlp609 <- read.csv(paste0(base_dir,'SA609_rna/deg_analysis/SA609-v6/SA609__fraction_dlp_corrected.csv'), check.names=F, stringsAsFactors=F)
  dlp1035 <- read.csv(paste0(base_dir,'SA1035_rna/deg_analysis/SA1035-v6/SA1035__fraction_dlp.csv'), check.names=F, stringsAsFactors=F)
  dlpSA535_cis <- read.csv(paste0(base_dir,'SA535_total_rna_v2/SA535-v6/SA535_SA535_cisplatin_fraction_dlp.csv'), check.names=F, stringsAsFactors=F)
  dlpSA535_cx <- read.csv(paste0(base_dir,'SA535_total_rna_v2/SA535-v6/SA535_SA535_CX5461_fraction_dlp.csv'), check.names=F, stringsAsFactors=F)
  
  # View(dlp609)
  # dlp609$desc_2 <- gsub('SA609_','',dlp609$desc)
  # dlp1035$desc_2 <- gsub('SA1035_','',dlp1035$desc)
  # dlpSA535_cis$desc_2 <- gsub('SA535_','',dlpSA535_cis$desc)
  # dlpSA535_cx$desc_2 <- gsub('SA535_','',dlpSA535_cx$desc)
  stat_dlp <- do.call(rbind,list(dlp609, dlp1035, dlpSA535_cis, dlpSA535_cx))
  # View(stat_dlp)
  stat_dlp$pct_incis_incr <- NULL
  stat_dlp$pct_incis_decr <- NULL
  stat_dlp$pct_others <- NULL 
  stat_dlp$pct_incis_incr <- round(as.numeric(stat_dlp$incis_incr)/as.numeric(stat_dlp$total_var_dlp) * 100, 2)
  stat_dlp$pct_incis_decr <- round(as.numeric(stat_dlp$incis_decr)/as.numeric(stat_dlp$total_var_dlp) * 100, 2)
  stat_dlp$pct_others <- 100 - stat_dlp$pct_incis_incr - stat_dlp$pct_incis_decr
  
  
  write.csv(stat_dlp, paste0(input_dir,'summary_dlp.csv'), row.names = F, quote = F)
  colnames(stat_dlp)
  stat_dlp1 <- stat_dlp[,c("desc","pct_incis_incr","pct_incis_decr","pct_others")]
  colnames(stat_dlp1) <- c("desc","Incis_incr","Incis_decr","Others_CNchange")
  
  stat_dlp1 <- stat_dlp1 %>%
    tidyr::pivot_longer(!desc, names_to = "gene_type", values_to = "percentage")
  
  dim(stat_dlp1)
  PDX <- lapply(strsplit(stat_dlp1$desc, "_"), function(x) {
    return(x[1])
  })
  unique(stat_dlp1$gene_type)
  stat_dlp1$PDX <- as.character(PDX)
  stat_dlp1$PDX <- ifelse(grepl('X',stat_dlp1$desc),'SA535_CX5461',
                          ifelse(!grepl('X',stat_dlp1$desc),paste0(stat_dlp1$PDX,'_CISPLATIN'),stat_dlp1$PDX))
  
  cols_dlp <- c("#E82503","#E61367","#00B300")
  stat_dlp1$gene_type <- factor(stat_dlp1$gene_type, levels = unique(stat_dlp1$gene_type))
  stat_dlp1$PDX <- factor(stat_dlp1$PDX, levels = unique(stat_dlp1$PDX))
  # stat_dlp1$desc <- gsub('SA609_','',stat_dlp1$desc)
  stat_dlp1$desc <- str_replace_all(stat_dlp1$desc, "(SA609_)|(SA535_)|(SA1035_)", "")
  p <- ggplot(stat_dlp1, aes(fill=gene_type, y=percentage, x=desc)) + 
    geom_bar(position="fill", stat="identity", width = 0.6) +
    scale_y_continuous(labels = scales::percent_format()) +
    facet_grid(. ~ PDX, scales="free", space='free') + 
    scale_fill_manual(values = cols_dlp) + 
    labs(x='DE Analysis: Resistant vs Sensitive cells', y="(%) genes in total CN change ", title='Proportion in-cis genes in total copy number change') +
    theme(legend.title = element_text(size=8), 
          legend.text = element_text(size=7),
          plot.title = element_text(color="black", size=13, hjust = 0.5),
          # legend.position = "none", 
          # axis.line = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          axis.text.x = element_text(size=8, angle = 90, color="black"),
          axis.text.y = element_text(size=8, color="black")
          #axis.title = element_blank(),
          # axis.ticks = element_blank()
    )
  
  png(paste0(output_dir,"cn_change_proportion.png"), height = 2*420, width=2*720, res = 2*72)
  print(p)
  dev.off()
  write.csv(stat_dlp,paste0(output_dir,'cn_change_proportion.csv'), row.names = F, quote = F)
  
  saveRDS(p, paste0(output_dir,"cn_change_proportion.rds"))
  
}

plot_intrans_incis_prevalence_v2 <- function(){
  
  # base_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/')
 
  gene_type=c('In_cis_Decrease_DownRegulated','In_cis_Decrease_UpRegulated',
              'In_cis_Increase_DownRegulated','In_cis_Increase_UpRegulated',
              'In_trans_DownRegulated','In_trans_UpRegulated'
  )
  gt <- data.frame(gene_type=gene_type,
                   Gene_Type=c('In cis Loss-Down','In cis Loss-Up',
                               'In cis Gain-Down','In cis Gain-Up',
                               'In trans Down','In trans Up'),
                   col_gt=c("#C3D7A4","#52854C",
                            "#FFDB6D","#D16103",
                            "#4E84C4","#293352"))
  
  keyvals_colour <- c("In cis Loss-Down"="#C3D7A4","In cis Loss-Up"="#52854C",
                      "In cis Gain-Down"="#FFDB6D","In cis Gain-Up"="#D16103",
                      "In trans Down"="#4E84C4","In trans Up"="#293352")
  
 
  stat <- stat %>% left_join(gt, by='gene_type')
  unique(stat$col_gt)
  # stat$DE_analysis <- gsub('(SA609_)|(SA1035_)|(SA535)','',stat$desc)
  # stat$DE_analysis <- gsub('^_','',stat$DE_analysis)
  write.csv(stat,paste0(save_dir,'summary_incis_intrans_genes_all_series.csv'), row.names = F, quote = F)
  # stat <- data.table::fread(paste0(save_dir,'summary_incis_intrans_genes_all_series.csv')) %>% as.data.frame()
  # input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/cis_trans_landscape_treated_and_untreated/differential_expression/'
  data_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/cis_trans/'
  pair_groups_fn <- paste0(data_dir,'comparisons_drug_res_v4.csv')
  pair_df <- data.table::fread(pair_groups_fn) %>% as.data.frame()
  
  rownames(pair_df) <- pair_df$desc
  pair_df <- pair_df %>%
    dplyr::select(-result_fn)
  
  summary(as.factor(stat$datatag))
  stat1 <- stat
  stat1 <- stat1 %>% left_join(pair_df, by=c('desc','datatag'))
  
  pdx_level = c('SA501','SA530','SA604','SA609','SA535','SA1035')
  stat1$PDX <- stat1$datatag
  stat1$PDX <- factor(stat1$PDX, levels = pdx_level)
  stat1$Gene_Type <- factor(stat1$Gene_Type, levels = gt$Gene_Type)
  # stat1$DE_desc1 <- gsub('M_N_O','MNO',stat1$DE_desc1)
  # stat1$DE_desc1 <- gsub('S_T','ST',stat1$DE_desc1)
  # summary(stat1$plt_desc)
  
  # t <- data.frame(plt_desc=unique(stat1$plt_desc))
  # data.table::fwrite(t, paste0(save_dir,'plt_desc.csv'), quote = F)
  # 
  my_font <- "Helvetica"
  # stat1$desc
  # stat1$plt_desc <- NULL
  save_dir <- input_dir
  data.table::fwrite(stat1, paste0(save_dir,'summary_incis_intrans_genes_plots.csv'), quote = F)
  
  
  # Plotting
  stat1 <- data.table::fread(paste0(save_dir,'summary_incis_intrans_genes_plots.csv')) %>% as.data.frame()
  plt_desc_df <- data.table::fread(paste0(save_dir,'plt_desc_v2.csv')) %>% as.data.frame()
  dim(plt_desc_df)
  dim(stat1)
  # View(plt_desc_df)
  
  plt_desc_df$ord <- rep(1:dim(plt_desc_df)[1],1)
  # stat1$plt_desc <- NULL
  # stat1$ord <- NULL
  # stat1 <- stat1 %>% inner_join(plt_desc_df, by=c('desc'))
  
  
  table(stat1$Gene_Type,stat1$PDX)
  rownames(stat1) <- paste0(stat1$desc, stat1$Gene_Type)
  p <- ggplot(stat1, aes(fill=Gene_Type, y=pct_genes, x=plt_desc)) + 
    geom_bar(stat="identity",position = "dodge",width = 0.8) +  #identity position="fill" , 
    scale_x_discrete(drop = FALSE) +
    # geom_text(aes(label=pct_genes), size=3.5) + #, color="white", vjust=1.6, 
    # geom_text(aes(label=ifelse(pct_genes < 5, NA, pct_genes)) , position = position_stack(vjust = 0.5), color="white", size=2.6) + #position=position_dodge(width=0.9), 
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    facet_grid(. ~ PDX, scales="free", space='free') + 
    scale_fill_manual(values = keyvals_colour, name='Gene Type')
  # p <- p + labs(x='DE Analysis: Resistant vs Sensitive cells', y="(%) Gene Type ", title='Differentially Expressed Genes - Copy Number Driven Proportion')
  p <- p + labs(x=NULL, y="(%) Gene Type ", title='Differentially Expressed Genes - Copy Number Driven Proportion')
  # p <- p + theme(text=element_text(family=my_font),
  #                legend.title = element_text(color="black", size=13, hjust = 0.5, family=my_font), 
  #                legend.text = element_text(color="black", size=11, hjust = 0.5, family=my_font),
  #                plot.title = element_text(color="black", size=14, hjust = 0.5, family=my_font, face = "bold"),
  #                legend.position = "bottom",
  #                # axis.line = element_blank(), 
  #                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
  #                panel.border = element_blank(),
  #                axis.text.x = element_text(size=12, angle = 90, color="black", family=my_font),
  #                axis.text.y = element_text(size=12, color="black", family=my_font),
  #                axis.title = element_text(size=12, color="black", family=my_font)
  #                # axis.ticks = element_blank()
  # )
  thesis_theme <- ggplot2::theme(
    text = element_text(size = 8, hjust = 0.5, family=my_font),
    axis.title.x = element_text(size=8, hjust = 0.5, family=my_font),  
    axis.title.y = element_text(size=8, hjust = 0.5, family=my_font),
    axis.text.x = element_text(size=7, hjust = 0.5, family=my_font, angle = 90),  
    axis.text.y = element_text(size=7, hjust = 0.5, family=my_font),
    plot.title = element_text(size=10, face="bold", hjust=0, family=my_font),
    strip.text.x = element_text(size=9),
    strip.text.y = element_text(size=9),
    legend.title=element_text(size=7, hjust = 0.5, family=my_font), 
    legend.text=element_text(size=7, hjust = 0.5, family=my_font),
    legend.spacing.x = unit(0.1, 'mm'),
    legend.spacing.y = unit(0.1, 'mm'),
    legend.key.height=unit(1,"line"),
    legend.position = "bottom",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
  )
  
  
  p <- p + thesis_theme
  output_dir <- save_dir
  png(paste0(output_dir,"incis_intrans_genes_prevalence_summary.png"), height = 2*600, width=2*650, res = 2*72)
  print(p)
  dev.off()
  ggsave(paste0(output_dir,"incis_intrans_genes_prevalence_summary.pdf"),
         plot = p,
         height = 5.5,
         width = 8,
         useDingbats=F)
  saveRDS(p, paste0(output_dir,"incis_intrans_genes_prevalence_summary.rds"))
  
  
  
  
}


plot_intrans_incis_prevalence_v4 <- function(stat, save_dir, verbose=F){

  stat <- stat %>%
    dplyr::filter(gene_type!='unmapped')
  
  stat$plt_desc <- paste0(stat$series,' ',stat$labels_detailed)
  stat$plt_desc <- gsub('vs.','vs.\n         ',stat$plt_desc)
  stat$desc <- stat$file_header
  my_font <- "Helvetica"
  plot_ls <- list()
  pd <- 'treated'
  df <- stat
  # df$gene_type
  # View(df)
  plot_ls[[pd]] <- dotplot_prevalence_v2(df, verbose = verbose, plttype='treated')
  # plot_ls[[pd]]#
  # dev.off()
  # df$desc <- df$order
  res_plt <- plot_cis_pos_neg(df, verbose = F, plttype='treated') # empty column for SA1035 here
  cistran_plt <- plot_cis_trans_barplot(df, verbose = verbose, plttype='treated')
  plot_ls[[paste0(pd,'_proportion')]] <- res_plt$p
  plot_ls[[paste0(pd,'_cistrans')]] <- cistran_plt$p
  ## Need 3 labels plot here, TO DO
  # View(head(tmp))
  for(tag in c('SA609','SA535','SA1035')){
    tmp <- df %>%
      dplyr::filter(datatag==tag)
    # colnames(tmp)
    plttitle1 <- NULL
    plttitle2 <- NULL
    plt_pg <- NULL
    # if(tag=='SA609'){
    #   plttitle1 <- 'Rx'
    #   plttitle2 <- 'UnRx'
    #   plt_pg <- 'Time'
    # }else{
    #   plttitle1 <- NULL
    #   plttitle2 <- NULL
    #   plt_pg <- NULL
    # }
    plot_ls[[paste0(pd,'_',tag,'_c1lb')]] <- viz_labels(tmp, lb_use='clone1',plottitle=plttitle1, #'Rx'
                                                color_scheme='predefined_clone_col',tag=tag)
    plot_ls[[paste0(pd,'_',tag,'_c2lb')]] <- viz_labels(tmp, lb_use='clone2',plottitle=plttitle2, #
                                                color_scheme='predefined_clone_col',tag=tag)
    plot_ls[[paste0(pd,'_',tag,'_passage')]] <- viz_labels(tmp, lb_use='passage',plottitle=plt_pg, 
                                                   color_scheme='gradient')
    
  }
  
  plot_ls[['lg']] <- get_legends_plt(df)
  plot_ls[['lg_cis']] <- res_plt$plg
  plot_ls[['lg_cistrans']] <- cistran_plt$plg
    
  return(plot_ls)
}  
plot_intrans_incis_prevalence_untreated_v4 <- function(stat, save_dir, verbose=F){
  stat <- stat %>%
    dplyr::filter(gene_type!='unmapped')
  
  stat$plt_desc <- paste0(stat$series,' ',stat$labels_detailed)
  stat$plt_desc <- gsub('vs.','vs.\n         ',stat$plt_desc)
  stat$desc <- stat$file_header
  my_font <- "Helvetica"
  plot_ls <- list()
  pd <- 'UnRx'
  df <- stat
  # df$gene_type
  # View(df)
  plot_ls[[pd]] <- dotplot_prevalence_v2(df, verbose = verbose, plttype=NULL)
  # plot_ls[[pd]]#
  # dev.off()
  # df$desc <- df$order
  res_plt <- plot_cis_pos_neg(df, verbose = verbose, plttype=NULL) # empty column for SA1035 here
  cistran_plt <- plot_cis_trans_barplot(df, verbose = verbose, plttype=NULL)
  plot_ls[[paste0(pd,'_proportion')]] <- res_plt$p
  plot_ls[[paste0(pd,'_cistrans')]] <- cistran_plt$p
  ## Need 3 labels plot here, TO DO
  plttitle1 <- NULL
  plttitle2 <- NULL
  plt_pg <- NULL
  pdxs_untreated = c('SA501','SA530','SA604','SA609','SA535','SA1035')
  for(tag in pdxs_untreated){
    tmp <- df %>%
      dplyr::filter(datatag==tag)
    plttitle1 <- NULL
    plttitle2 <- NULL
    plt_pg <- NULL
    plot_ls[[paste0(pd,'_',tag,'_c1lb')]] <- viz_labels(tmp, lb_use='clone1',plottitle=plttitle1, #'Rx'
                                                        color_scheme='predefined_clone_col',tag=tag)
    plot_ls[[paste0(pd,'_',tag,'_c2lb')]] <- viz_labels(tmp, lb_use='clone2',plottitle=plttitle2, #
                                                        color_scheme='predefined_clone_col',tag=tag)
    plot_ls[[paste0(pd,'_',tag,'_passage')]] <- viz_labels(tmp, lb_use='passage',plottitle=plt_pg, 
                                                           color_scheme='gradient',tag=tag)
  }
    
  plot_ls[['lg']] <- get_legends_plt(df)
  plot_ls[['lg_cis']] <- res_plt$plg
  plot_ls[['lg_cistrans']] <- cistran_plt$plg
  
  return(plot_ls)
}  

plot_intrans_incis_prevalence_v3 <- function(stat, verbose=F){
  
  # meta_ptx <- data.frame(datatag=c("SA501","SA530","SA604","SA609","SA535","SA1035"), 
  #                        pt=paste0("Pt",rep(1:6,1)))
  # stat <- stat %>% left_join(meta_ptx, by='datatag')
  # dim(stat)

  # stat$plt_desc <- paste0(stat$pt,':\n',stat$plt_desc)
  stat$plt_desc <- paste0(stat$series,' ',stat$labels_detailed)
  stat$plt_desc <- gsub('vs.','vs.\n         ',stat$plt_desc)
  stat$desc <- stat$file_header
  pdxs = c('SA609','SA535','SA1035')
  pdxs_untreated = c('SA501','SA530','SA604','SA609','SA535','SA1035')
  my_font <- "Helvetica"
  
  # First get untreated comparison
  df_untreated <- stat %>% 
    dplyr::filter(datatag %in% pdxs_untreated & comp_type=='untreated')
  
  plot_ls <- list()
  # p <- dotplot_prevalence_v2(df_untreated, verbose = verbose)
  # p <- dotplot_prevalence(df_untreated, verbose = verbose)
  # df_untreated$plt_desc <- factor(df_untreated$plt_desc, levels = unique(df_untreated$plt_desc))
  # Untreated series ## TO DO: add patient names here
  plot_ls[['untreated']] <- dotplot_prevalence_v2(df_untreated, verbose = verbose)
  # df_untreated$desc <- df_untreated$file_header
  df_untreated$desc <- df_untreated$order
  res_untreated_plt <- plot_cis_pos_neg(df_untreated, verbose = verbose)
  res_untreated_cistran_plt <- plot_cis_trans_barplot(df_untreated, verbose = verbose)
  plot_ls[['untreated_proportion']] <- res_untreated_plt$p
  plot_ls[['untreated_cistrans']] <- res_untreated_cistran_plt$p
  plot_ls[['untreated_c1lb']] <- viz_labels(df_untreated, lb_use='clone1',plottitle='UnRx', color_scheme='untreated',tag='untreated')
  plot_ls[['untreated_c2lb']] <- viz_labels(df_untreated, lb_use='clone2',plottitle='UnRx', color_scheme='untreated',tag='untreated')
  
  df_untreated$passage <- 'X' ## TO DO
  plot_ls[['untreated_passage']] <- viz_labels(df_untreated, lb_use='passage',plottitle='Time', color_scheme='gradient')
  
  # Second, for time series
  for(pd in pdxs){
    df <- stat %>% 
      dplyr::filter(datatag==pd & comp_type=='treated_vs_untreated')
    plot_ls[[pd]] <- dotplot_prevalence_v2(df, verbose = verbose)
    df$desc <- df$order
    res_plt <- plot_cis_pos_neg(df, verbose = verbose)
    cistran_plt <- plot_cis_trans_barplot(df, verbose = verbose)
    plot_ls[[paste0(pd,'_proportion')]] <- res_plt$p
    plot_ls[[paste0(pd,'_cistrans')]] <- cistran_plt$p
    plot_ls[[paste0(pd,'_c1lb')]] <- viz_labels(df, lb_use='clone1',plottitle='Rx', color_scheme='predefined_clone_col',tag='SA609')
    plot_ls[[paste0(pd,'_c2lb')]] <- viz_labels(df, lb_use='clone2',plottitle='UnRx', color_scheme='predefined_clone_col',tag='SA609')
    if(pd=='SA535'){ ## issue with meta file, TO DO
      df$passage <- ifelse(df$passage=='X1','X10',df$passage)
    }
    plot_ls[[paste0(pd,'_passage')]] <- viz_labels(df, lb_use='passage',plottitle='Time', color_scheme='gradient')
    
    
  }
  
  plot_ls[['lg']] <- get_legends_plt(df)
  plot_ls[['lg_cis']] <- res_untreated_plt$plg
  plot_ls[['lg_cistrans']] <- res_untreated_cistran_plt$plg
  return(plot_ls)
}  
get_legends_plt <- function(df){
  # keyvals_colour <- c("In cis Loss-Down"="#C3D7A4","In cis Loss-Up"="#52854C",
  #                     "In cis Gain-Down"="#FFDB6D","In cis Gain-Up"="#D16103",
  #                     "In trans Down"="#4E84C4","In trans Up"="#293352")
  keyvals_colour <- c("Loss-Down "="#C3D7A4","Loss-Up "="#52854C",
                      "Gain-Down "="#FFDB6D","Gain-Up "="#D16103",
                      "Trans Down "="#4E84C4","Trans Up "="#293352")
  # gt_displ <- gsub('In cis ','',names(keyvals_colour))
  # gt_displ <- gsub('In t','T',gt_displ)
  gt <- c("In_cis_Decrease_DownRegulated","In_cis_Decrease_UpRegulated",
          "In_cis_Increase_DownRegulated","In_cis_Increase_UpRegulated",
          "In_trans_DownRegulated","In_trans_UpRegulated")
  
  gt1 <- names(keyvals_colour)
  names(gt1) <- gt
  df$Gene_Type <- gt1[df$gene_type]
  p <- ggplot(df, aes(x = pct_genes, y = plt_desc, color = Gene_Type)) + 
    geom_point(size=5, alpha=1) + #, size=2*log10(pct_genes)
    scale_color_manual(name = NULL, values = keyvals_colour) + 
    theme(legend.position ='bottom',
          legend.text = element_text(size = 12, hjust = 0.5, family=my_font, angle = 90),
          legend.title = element_text(size = 10, hjust = 0.5, family=my_font, angle = 90))
  p <- p + guides(color = guide_legend(override.aes = list(size=8, ncol=1)))#shape = 0, , ncol = 2
  lg <- cowplot::get_legend(p)
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  return(plg)
}  

dotplot_prevalence_v2 <- function(df, verbose=F, plttype = NULL){
  my_font <- "Helvetica"
  keyvals_colour <- c("In cis Gain-Up"="#D16103","In cis Loss-Down"="#C3D7A4",
                      "In cis Gain-Down"="#FFDB6D","In cis Loss-Up"="#52854C",
                      "In trans Down"="#4E84C4","In trans Up"="#293352")
  keyvals_colour <- c("Cis Gain-Up"="#D16103","Cis Loss-Down"="#C3D7A4",
                      "Cis Gain-Down"="#FFDB6D","Cis Loss-Up"="#52854C",
                      "Trans Down"="#4E84C4","Trans Up"="#293352")
  
  gt <- c("In_cis_Increase_UpRegulated","In_cis_Decrease_DownRegulated",
          "In_cis_Increase_DownRegulated","In_cis_Decrease_UpRegulated",
          "In_trans_DownRegulated","In_trans_UpRegulated")
  
  
  gt1 <- names(keyvals_colour)
  names(gt1) <- gt
  df$Gene_Type <- gt1[df$gene_type]
  # unique(df$Gene_Type)
  df <- df %>%
    dplyr::arrange(-order)
  # df1 <- df[!duplicated(df$plt_desc),]
  df$plt_desc <- factor(df$plt_desc, levels = unique(df$plt_desc))
  df$Gene_Type <- factor(df$Gene_Type, levels = names(keyvals_colour))
  # df$sz <- ifelse(df$pct_genes>0, log2(df$pct_genes)+1, 0.01)
  df$sz <- ifelse(log2(df$pct_genes)+1>0, log2(df$pct_genes)+1, 0.01)
  # max_pct <- log2(max(df$pct_genes))+1
  # print(max_pct)
  p <- ggplot(df, aes(x=Gene_Type, y=plt_desc)) + 
    geom_point(aes(color = Gene_Type, size=sz), alpha=1) + #, size=2*log10(pct_genes), size=4.5
    scale_color_manual(values = keyvals_colour, guide = "none") +
    guides(size=guide_legend(title="Log2(%gene)+1 ")) #+ 
    # lims(size=c(NA,max_pct))
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    # scale_x_continuous(trans = scales::log2_trans(),
    #                    breaks = scales::trans_breaks("log2", function(x) 2^x))+ #, labels=scales::percent_format(scale = 1)
    # breaks = scales::pretty_breaks(n = 5), labels=scales::percent
    # scale_y_discrete(position = "right")
  # scale_x_continuous(breaks = scales::log2_trans()) + 
  # facet_grid(rows = vars(PDX)) + 
  # geom_hline(yintercept = 1:6, col = "#e2e2e2") +
  # geom_text(aes(label=celltype_desc))+
  # annotate('text', x = df$success, y = df$index, label = df$success, size=3, colour = col[df$gender])+
  if(!is.null(plttype)){
    p <- p + facet_grid(series ~ ., scales="free_y", space='free')
  }
  # scale_color_manual(values = col) +
  # p1 <- p
  # p <- p + guides(fill = guide_legend(override.aes = list(shape = 0, size=6, nrow = 2)))
  # lg <- cowplot::get_legend(p)
  # plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  
  p <- p + labs(title=NULL, x=NULL, y=NULL) + #'Gene type'
    ct_theme
    # theme_bw(base_size = 10) + 
    # theme(#legend.position = "bottom",
    #       legend.position = "none",
    #       # panel.grid.major.x = element_blank(),
    #       # panel.grid.minor.x = element_blank(),
    #       # axis.ticks.y = element_blank(),
    #       # axis.text  = element_blank(),
    #       text = element_text(size = 11, hjust = 0.5, family=my_font),
    #       # axis.text.x = element_text(size=11, hjust = 0.5, family=my_font, angle = 90),  #
    #       axis.text.x = element_blank(),
    #       # axis.text.x = element_text(size=11, hjust = 0.5, family=my_font, angle = 90, color='black'), ## for testing
    #       # axis.text.y = element_text(size=11, hjust = 0.5, family=my_font), ## for testing
    #       axis.text.y = element_blank(),
    #       plot.title = element_text(size=11, hjust=0.5, family=my_font), #, face="bold"
    #       axis.title.y = element_blank(),
    #       # axis.title.x = element_text(size=11, hjust = 0.5, family=my_font)
    #       # axis.title.x = element_blank()
    #       # strip.text.x = element_text(color="black",size=11, family=my_font),
    #       # strip.text.y = element_text(color="black",size=11, family=my_font),
    #       strip.text = element_text(size=11, color="black", family=my_font),
    #       strip.background = element_blank(),
    #       # strip.placement = "outside",
    #       panel.background = element_rect(fill = "white", colour = "grey50")
    #     ) #+
  # p
  # xlim(0,xmax)
  # geom_text_repel( data = df, aes(label = pct_genes), nudge_x=0, nudge_y=0, max.overlaps = Inf,
  #                  min.segment.length = 0,segment.size=0.1,
  #                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines") )
  return(p)
}

dotplot_prevalence <- function(df, verbose=F){
  my_font <- "Helvetica"
  keyvals_colour <- c("In cis Loss-Down"="#C3D7A4","In cis Loss-Up"="#52854C",
                      "In cis Gain-Down"="#FFDB6D","In cis Gain-Up"="#D16103",
                      "In trans Down"="#4E84C4","In trans Up"="#293352")
  gt <- c("In_cis_Decrease_DownRegulated","In_cis_Decrease_UpRegulated",
                      "In_cis_Increase_DownRegulated","In_cis_Increase_UpRegulated",
                      "In_trans_DownRegulated","In_trans_UpRegulated")
  
  gt1 <- names(keyvals_colour)
  names(gt1) <- gt
  df$Gene_Type <- gt1[df$gene_type]
  # unique(df$Gene_Type)
  df <- df %>%
    dplyr::arrange(-order)
  # df1 <- df[!duplicated(df$plt_desc),]
  df$plt_desc <- factor(df$plt_desc, levels = unique(df$plt_desc))
  df$sz <- ifelse(df$pct_genes>0, 2*log10(df$pct_genes), 0.01)
  p <- ggplot(df, aes(x = pct_genes, y = plt_desc)) + 
    geom_point(aes(color = Gene_Type, size=sz), alpha=0.9) + #, size=2*log10(pct_genes), size=4.5
    scale_color_manual(values = keyvals_colour) + 
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    scale_x_continuous(trans = scales::log2_trans(),
                       breaks = scales::trans_breaks("log2", function(x) 2^x))+ #, labels=scales::percent_format(scale = 1)
    # breaks = scales::pretty_breaks(n = 5), labels=scales::percent
    scale_y_discrete(position = "right")
    # scale_x_continuous(breaks = scales::log2_trans()) + 
    # facet_grid(rows = vars(PDX)) + 
    # geom_hline(yintercept = 1:6, col = "#e2e2e2") +
    # geom_text(aes(label=celltype_desc))+
    # annotate('text', x = df$success, y = df$index, label = df$success, size=3, colour = col[df$gender])+
    
    # scale_color_manual(values = col) +
  p <- p + labs(title=' ', x='Gene type prevalence (%)') +
    theme_bw(base_size = 10) + 
    theme(legend.position = "none",
          # panel.grid.major.x = element_blank(),
          # panel.grid.minor.x = element_blank(),
          # axis.ticks.y = element_blank(),
          # axis.text  = element_blank(),
          text = element_text(size = 8, hjust = 0.5, family=my_font),
          axis.text.x = element_text(size=9, hjust = 0.5, family=my_font),  #, angle = 90
          axis.text.y = element_text(size=9, hjust = 0.5, family=my_font),
          plot.title = element_text(size=10, hjust=0.5, family=my_font), #, face="bold"
          axis.title.y = element_blank(),
          axis.title.x = element_text(size=9, hjust = 0.5, family=my_font)) #+
  # p
  # xlim(0,xmax)
  # geom_text_repel( data = df, aes(label = pct_genes), nudge_x=0, nudge_y=0, max.overlaps = Inf,
  #                  min.segment.length = 0,segment.size=0.1,
  #                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines") )
  return(p)
}

plot_cis_trans_barplot <- function(df, verbose=F, plttype=NULL){
  my_font <- "Helvetica"
  GT <- c("#D2691E","blue2") #chocolate -cis
  names(GT) <- c("In cis","In trans") #
  stat1 <- df %>%
    dplyr::mutate(gt=case_when(
      grepl('In_cis',gene_type) ~ "In cis",
      grepl('In_trans',gene_type) ~ "In trans",
      TRUE ~ "unmapped"
    )) %>%
    dplyr::group_by(gt, order,series)%>% 
    dplyr::summarise(pct_genes=sum(pct_genes))#%>% 
    # dplyr::select(gt,pct_genes,desc)%>% 
    # tidyr::pivot_wider(names_from = 'gt', values_from = 'pct_genes',values_fill=0)# # df[is.na(df)] <- 0
  # dim(stat1)
  
  # stat1$desc <- factor(stat1$desc, levels = unique(stat1$desc))
  stat1$gt <- factor(stat1$gt, levels = c("In trans","In cis"))
  # summary(stat3$pct_genes)
  p <- ggplot(stat1, aes(fill=gt, y=pct_genes, x=order)) + 
    geom_bar(position="fill", stat="identity", width = 0.5) + 
    coord_flip() + 
    scale_fill_manual(values = GT, name='',labels = names(GT)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_x_discrete(drop=T)
  
  
  p <- p +
    # ct_theme
     theme(text=element_text(family=my_font),
          axis.text.x = element_text(size=11, color="black", family=my_font),
          axis.text.y = element_text(size=9, color="black", family=my_font),
          axis.title.y = element_blank(),
          # legend.position ='bottom',
          plot.title = element_text(size=10, hjust=0.5, family=my_font),
          legend.text = element_text(size=9, color="black", family=my_font),
          legend.title = element_blank())
  p <- p + guides(fill = guide_legend(override.aes = list(shape = 0, size=6, nrow = 2)))
  lg <- cowplot::get_legend(p)
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  
  if(verbose==F){
    axisY <- ggplot2::element_blank()
  }else{
    axisY <- ggplot2::element_text(size=9, hjust = 0.5, family=my_font)
  }
  p1 <- ggplot(stat1, aes(fill=gt, y=pct_genes, x=order)) + 
    geom_bar(position="fill", stat="identity", width = 0.7) + 
    coord_flip() + 
    scale_fill_manual(values = GT, name='Gene Type',labels = names(GT)) + #, drop=F
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +#, labels=scales::percent
    scale_x_discrete(drop=T)
  if(!is.null(plttype)){
    p1 <- p1 + facet_grid(series ~ ., scales="free_y", space='free')
  }
  p1 <- p1 + labs(title=NULL, y=NULL, x=NULL) + #'Cis/trans genes fraction'
    ct_theme
    # theme(plot.title = element_text(size=11, hjust=0.5, family=my_font),
    #       text=element_text(family=my_font),
    #       legend.position = 'none',
    #       axis.text.x = element_text(size=11, hjust = 0.5, family=my_font),
    #       axis.title.x = element_text(size=11, hjust = 0.5, family=my_font),
    #       # axis.text.y = element_blank(),
    #       axis.text.y = axisY,
    #       axis.title.y = element_blank(),
    #       # axis.ticks.length=unit(0.25, "cm"),
    #       # axis.ticks.y = element_blank(),
    #       strip.text.y = element_blank(), 
    #       panel.background = element_rect(fill = 'white', colour = 'white'))
  
  return(list(plg=plg, p=p1))
}
plot_cis_pos_neg <- function(df, verbose=F, plttype=NULL){
  my_font <- "Helvetica"
  cis_cols <- c('lightseagreen','red')
  # cis_cols <- c("chocolate","light green")
  
  names(cis_cols) <- c("incis_positive","incis_negative")
  lbcis <- c("In cis canonical", 
             "In cis non canonical")
  
  # excluded <- df %>% 
  #   dplyr::mutate(exist_cis=grepl('In_cis',gene_type)) %>% 
  #   dplyr::group_by(file_header, order, series) %>% 
  #   dplyr::summarise(nb_cis=sum(exist_cis)) %>% 
  #   dplyr::filter(nb_cis==0) %>% 
  #   as.data.frame()
    # dplyr::pull(file_header)
  
  df <- df %>% 
    dplyr::filter(grepl('In_cis',gene_type))%>% 
    dplyr::select(gene_type,pct_genes,order, series)%>% 
    tidyr::pivot_wider(names_from = 'gene_type', values_from = 'pct_genes',values_fill=0)# # df[is.na(df)] <- 0

  df <- df %>%
    dplyr::mutate(incis_positive=In_cis_Decrease_DownRegulated+In_cis_Increase_UpRegulated,
                  incis_negative=In_cis_Decrease_UpRegulated+In_cis_Increase_DownRegulated)
  stat3 <- df %>%
    dplyr::select(incis_positive, incis_negative, order, series)
  # tibble::column_to_rownames(desc)%>% 
  # stat3 <- stat3 %>% tibble::column_to_rownames('desc')
  # head(stat3)
  stat3 <- stat3 %>% tidyr::pivot_longer(cols = c('incis_positive','incis_negative'), 
                                         names_to = "gene_type", values_to = "pct_genes") %>%
    dplyr::arrange(-order)
  # stat3$order <- factor(stat3$order, levels = unique(stat3$order))
  # summary(stat3$pct_genes)
  # View(stat3)
  # colnames(stat3)
  # rownames(excluded) <- excluded$file_header
  # for(fh in excluded$file_header){
  #   tmp <- data.frame(order=rep(excluded[fh,'order'],2), series=rep(excluded[fh,'series'],2),
  #                     gene_type=c('incis_positive','incis_negative'), pct_genes=c(0,0))
  #   stat3 <- dplyr::bind_rows(stat3, tmp)
  # }
  
  p <- ggplot(stat3, aes(fill=gene_type, y=pct_genes, x=order)) + 
    geom_bar(position="fill", stat="identity", width = 0.5) + 
    coord_flip() + 
    scale_fill_manual(values = cis_cols, name='',labels = lbcis) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) #+ 
    
  p <- p + theme(text=element_text(family=my_font),
          axis.text.x = element_text(size=9, color="black", family=my_font),
          axis.text.y = element_text(size=9, color="black", family=my_font),
          axis.title.y = element_blank(),
          # legend.position ='bottom',
          plot.title = element_text(size=10, hjust=0.5, family=my_font),
          legend.text = element_text(size=9, color="black", family=my_font),
          legend.title = element_blank())
  p <- p + guides(fill = guide_legend(override.aes = list(shape = 0, size=6, nrow = 2)))
  lg <- cowplot::get_legend(p)
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  
  if(verbose==F){
    axisY <- ggplot2::element_blank()
  }else{
    axisY <- ggplot2::element_text(size=9, hjust = 0.5, family=my_font)
  }
  # stat3$pct_genes <- stat3$pct_genes * 100
  # View(stat3)
  # stat3$order
  p1 <- ggplot(stat3, aes(fill=gene_type, y=pct_genes, x=order)) + 
    geom_bar(position="fill", stat="identity", width = 0.7) + 
    coord_flip() + 
    scale_fill_manual(values = cis_cols, name='Gene Type',labels = lbcis) +#+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +#, labels=scales::percent
    scale_x_discrete(drop=T)
  # p1
  if(!is.null(plttype)){
    p1 <- p1 + facet_grid(series ~ ., scales="free_y", space='free')
  }
  p1 <- p1 + labs(title=NULL, y=NULL, x=NULL) + #'Fraction cis gene', x=NULL 'Cis canonical/non canonical'
    ct_theme
    # theme(plot.title = element_text(size=11, hjust=0.5, family=my_font),
    #       text=element_text(family=my_font),
    #       legend.position = 'none',
    #       axis.text.x = element_text(size=11, hjust = 0.5, family=my_font),
    #       axis.title.x = element_text(size=11, hjust = 0.5, family=my_font),
    #       # axis.text.y = element_blank(),
    #       axis.text.y = axisY,
    #       axis.title.y = element_blank(),
    #       # axis.ticks.length=unit(0.25, "cm"),
    #       # axis.ticks.y = element_blank(),
    #       strip.background = element_blank(),
    #       strip.text.y = element_blank(),
    #       panel.background = element_rect(fill = 'white', colour = 'white')
    #       )
  
  # p1 <- p1 + labs(title='Cis/trans genes fraction', y=NULL) + 
  #   theme(plot.title = element_text(size=11, hjust=0.5, family=my_font),
  #         text=element_text(family=my_font),
  #         legend.position = 'none',
  #         axis.text.x = element_text(size=11, hjust = 0.5, family=my_font),
  #         axis.title.x = element_text(size=11, hjust = 0.5, family=my_font),
  #         # axis.text.y = element_blank(),
  #         axis.text.y = axisY,
  #         axis.title.y = element_blank(),
  #         # axis.ticks.length=unit(0.25, "cm"),
  #         # axis.ticks.y = element_blank(),
  #         strip.text.y = element_blank(), 
  #         panel.background = element_rect(fill = 'white', colour = 'white'))
  # p1
  return(list(plg=plg, p=p1))
}

plot_all_volcano <- function(base_dir, input_dir){
  input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results','/SA609_rna/deg_analysis/SA609-v6/')
  p609 <- readRDS(paste0(input_dir,'SA609_UTTT_R_UUUU_C/DE_SA609_UTTT_R_UUUU_C.rds'))
  input_dir <- paste0(base_dir,'SA1035_rna/deg_analysis/SA1035-v6/')
  p1035 <- readRDS(paste0(input_dir,'SA1035_UTTTT_G_UUUUU_E/DE_SA1035_UTTTT_G_UUUUU_E.rds'))
  input_dir <- paste0(base_dir,'SA535_total_rna_v2/SA535-v6/')
  pSA535_cis <- readRDS(paste0(input_dir,'SA535_UUTTT_S_T_UUUUU_Q/DE_SA535_UUTTT_S_T_UUUUU_Q.rds'))
  pSA535_cx <- readRDS(paste0(input_dir,'SA535_UXXXX_U_UUUUU_Q/DE_SA535_UXXXX_U_UUUUU_Q.rds'))
  save_dir <- input_dir
  main_plot <- cowplot::plot_grid(p609, p1035, pSA535_cis, pSA535_cx,
                                  ncol = 2,
                                  nrow=2,
                                  # rel_heights = c(2.5,1),
                                  align = 'hv'
  )
  png(paste0(save_dir,"volcano_3series.png"), height = 2*920, width=2*1050, res = 2*72)
  print(main_plot)
  dev.off()
}

plot_fitness_genes <- function(){
  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
  input_dir <- paste0(base_dir,'SA535_total_rna_v2/SA535-v6/')
  pair_groups_fn <- paste0(input_dir,'comparisons_drug_res.csv')
  id609 <- paste0(base_dir,'SA609_rna/deg_analysis/SA609-v6/')
  
  id1035 <- paste0(base_dir,'SA1035_rna/deg_analysis/SA1035-v6/')
  id535 <- paste0(base_dir,'SA535_total_rna_v2/SA535-v6/')
  pair_groups <- read.csv(pair_groups_fn, header=T, check.names=F, stringsAsFactors=F)
  View(pair_groups)
  dim(pair_groups)
  length(unique(pair_groups$desc))
  unique(pair_groups$datatag)
  pair_groups$result_dir <- ifelse(pair_groups$datatag=='SA609',paste0(id609,pair_groups$desc,'/'),
                                   ifelse(pair_groups$datatag=="SA1035",paste0(id1035,pair_groups$desc,'/'),
                                          paste0(id535,pair_groups$desc,'/')))
  
  if(is.null(cancer_ref_genes_fn)){
    ref_dif <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/'
    cancer_ref_genes_fn <- paste0(ref_dif,'Behan_CFgenes.csv')
  }
  cancer_ref_genes_df <- read.csv(cancer_ref_genes_fn, stringsAsFactors=F, check.names = F)
  print(dim(cancer_ref_genes_df))
  colnames(cancer_ref_genes_df)[which(colnames(cancer_ref_genes_df) == "ADAM PanCancer Core-Fitness genes")] <- "PanCancer_Fitness_genes"
  
  
  cancer_genes_ls <- list()
  rownames(pair_groups) <- pair_groups$desc
  minLogFC <- 0.25
  FDR_cutoff <- 0.01
  pValueThrs <- 0.05
  for(de in pair_groups$desc){
    signif_genes <- read.csv(paste0(pair_groups[de,'result_dir'],'signif_genes.csv'), check.names = F, stringsAsFactors=F)
    signif_genes <- signif_genes[abs(signif_genes$logFC)>minLogFC & signif_genes$FDR<FDR_cutoff  & signif_genes$PValue<pValueThrs,]
    signif_genes <- signif_genes[signif_genes$is_fitness_gene,]
    print(dim(signif_genes))
    # fitness_genes <- summary(as.factor(signif_genes$classified_gene_dlp))
    nb_incis <- sum(!is.na(signif_genes$classified_gene_dlp)==TRUE)
    nb_intrans <- sum(is.na(signif_genes$classified_gene_dlp)==TRUE)
    # cancer_genes_ls[[paste0(datatag,'_',de)]] <- round(as.numeric(fitness_genes['TRUE'])*100/nrow(cancer_ref_genes_df),2)
    cancer_genes_ls[[paste0(de,'_incis')]] <- round(nb_incis*100/nrow(cancer_ref_genes_df),2)
    cancer_genes_ls[[paste0(de,'_intrans')]] <- round(nb_intrans*100/nrow(cancer_ref_genes_df),2)
  }
  summary(as.factor(pair_groups$title))
  unique(pair_groups$title)
  cancer_pct_df <- do.call(rbind, cancer_genes_ls)
  cancer_pct_df <- as.data.frame(cancer_pct_df)
  dim(cancer_pct_df)
  View(head(cancer_pct_df))
  colnames(cancer_pct_df)[which(colnames(cancer_pct_df)=='V1')] <- 'Fitness_genes_pct'
  series = lapply(strsplit(rownames(cancer_pct_df), "_"), function(x) {
    return(x[1])
  })
  cancer_pct_df$PDX <- as.character(series)
  
  level_pdx <- c('SA609_CISPLATIN','SA1035_CISPLATIN','SA535_CISPLATIN','SA535_CX5461')
  cancer_pct_df$PDX <- ifelse(cancer_pct_df$PDX=='SA609','SA609_CISPLATIN',
                              ifelse(cancer_pct_df$PDX=='SA1035','SA1035_CISPLATIN',cancer_pct_df$PDX))
  cancer_pct_df$PDX <- ifelse(grepl('SA535',rownames(cancer_pct_df)) & grepl('T',rownames(cancer_pct_df)),'SA535_CISPLATIN',
                              ifelse(grepl('SA535',rownames(cancer_pct_df)) & grepl('X',rownames(cancer_pct_df)),'SA535_CX5461',cancer_pct_df$PDX))
  summary(as.factor(cancer_pct_df$PDX))
  View(rownames(cancer_pct_df))
  # cancer_pct_df$PDX <- c(rep(level_pdx, c(12, 8, 14, 12)))
  cancer_pct_df$PDX <- factor(cancer_pct_df$PDX, levels=level_pdx)
  
  
  rn <- as.character(rownames(cancer_pct_df))
  cancer_pct_df$desc <- rn #factor(rn, levels = rn)
  output_dir <- input_dir
  write.csv(cancer_pct_df, paste0(output_dir,'total_cancer_pct_df.csv'), quote=F, row.names = F)
  cancer_pct_df <- read.csv(paste0(output_dir,'total_cancer_pct_df.csv'), check.names = F, stringsAsFactors = F)
  View(cancer_pct_df)
  cancer_pct_df$gene_type <- ifelse(grepl("*incis",cancer_pct_df$desc),"incis",'intrans')
  # viz_pancancer_genes(cancer_pct_df, output_dir, tag="totaldata")
  unique(cancer_pct_df$gene_type)
  # desc <- mclapply(strsplit(cancer_pct_df$desc, "_"), function(x) {
  #   if(length(x)==7){
  #     return(paste(x[2],x[3],x[4],x[5], sep='_'))
  #   }else{
  #     return(paste(x[2],x[3],x[4],x[5],x[6], sep='_'))
  #   }
  #   
  # }, mc.cores = 2)
  cancer_pct_df$de_pair <- cancer_pct_df$desc
  # cancer_pct_df$de_pair <- gsub('_incis', '', grep("*incis",cancer_pct_df$de_pair, value=T))
  # cancer_pct_df$de_pair <- gsub('_intrans', '', grep("*intrans",cancer_pct_df$de_pair, value=T))
  cancer_pct_df$de_pair <- ifelse(grepl("*incis",cancer_pct_df$de_pair), gsub('_incis', '', cancer_pct_df$de_pair),
                                  ifelse(grepl("*intrans",cancer_pct_df$de_pair), gsub('_intrans', '', cancer_pct_df$de_pair),cancer_pct_df$de_pair))
  # View(cancer_pct_df$de_pair)
  cancer_pct_df$de_pair <- gsub('^(SA609_)|(SA1035_)|(SA535_)','',cancer_pct_df$de_pair)
  colorcode <- c('#FF3232','#40A0E0')
  ylabel <- '(%) Core-Fitness genes '
  xlabel <- 'DE Analysis: Resistant versus Sensitive cells'
  plottitle <- 'Percentage genes in ADAM PanCancer Core-Fitness genes'
  
  # cancer_pct_df$de_pair <- factor(cancer_pct_df$de_pair, levels = unique(cancer_pct_df$de_pair))
  
  # cancer_pct_df$series <- factor(cancer_pct_df$series, levels = unique(cancer_pct_df$series))
  plot_stack_barplot(cancer_pct_df, colorcode, xlabel, ylabel, plottitle, output_dir,'totaldata',
                     fa='gene_type', xa='de_pair', ya='Fitness_genes_pct')
  
  
  
  # Genes that appear in 20 cisplatin genes in 3 series
  rownames(pair_groups) <- pair_groups$desc
  minLogFC <- 0.25
  FDR_cutoff <- 0.01
  pValueThrs <- 0.05
  de_genes_ls <- list()
  for(de in pair_groups$desc){
    signif_genes <- read.csv(paste0(pair_groups[de,'result_dir'],'signif_genes.csv'), check.names = F, stringsAsFactors=F)
    signif_genes <- signif_genes[abs(signif_genes$logFC)>minLogFC & signif_genes$FDR<FDR_cutoff  & signif_genes$PValue<pValueThrs,]
    # signif_genes <- signif_genes[signif_genes$is_fitness_gene,]
    print(dim(signif_genes))
    colnames(signif_genes)
    signif_genes <- signif_genes[,c("gene_symbol","ensembl_gene_id")]
    signif_genes$datatag <- pair_groups[de,'datatag']
    dim(signif_genes)
    de_genes_ls[[paste0(de)]] <- signif_genes
  }
  summary(as.factor(pair_groups$title))
  unique(pair_groups$title)
  de_genes_df <- do.call(rbind, de_genes_ls)
  de_genes_df <- as.data.frame(de_genes_df)
  View(head(de_genes_df))
  dim(de_genes_df)
  data.table::fwrite(de_genes_df,file = paste0(input_dir,'total_DE_genes_3series.csv'), sep=',', quote=F)
  de_genes_df$description <- paste(de_genes_df$gene_symbol,de_genes_df$datatag, sep='_')
  de_genes_df1 <- de_genes_df[!duplicated(de_genes_df$description),]
  dim(de_genes_df1)
  View(head(de_genes_df1))
  colnames(de_genes_df1)
  de_genes_df2 <- de_genes_df1[,c("gene_symbol","datatag")]
  # de_genes_df2 <- de_genes_df2 %>%
  #   pivot_wider(names_from = datatag, values_from = gene_symbol)
  
  
  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
  cis_20 <- read.csv(paste0(base_dir,'biodatabase/cisplatin_20_genes.csv'), check.names = F, stringsAsFactors=F)
  
  dim(cis_20)
  de_genes_20 <- de_genes_df2[de_genes_df2$gene_symbol %in% cis_20$gene,]
  dim(de_genes_20)
  
  rownames(cis_20) <- cis_20$gene
  cis_20$DE_expression_series <- NA
  for(g in cis_20$gene){
    c=de_genes_20[de_genes_20$gene_symbol==g,'datatag']
    print(paste(c,collapse = '_'))
    if(!is.null(c)){
      cis_20[g,'DE_expression_series'] <- paste(c,collapse = '_')
    }
  }
  
  write.csv(cis_20, file = paste0(input_dir,'cisplatin_20genes_3series.csv'), quote=F)
  colnames(de_genes_20)
  # rownames(de_genes_20) <- de_genes_20
  de_genes_20$gene_symbol <- as.character(de_genes_20$gene_symbol)
  de_genes_20 <- de_genes_20 %>%
    pivot_wider(names_from = gene_symbol, values_from = datatag, values_fn = list)
  View(de_genes_20)
  de_genes_20 <- t(de_genes_20)
  de_genes_20 <- as.data.frame(de_genes_20)
  colnames(de_genes_20) <- 'DE_expression_series'
  de_genes_20$gene_symbol <- rownames(de_genes_20)
  de_genes_20$DE_expression_series <- paste(de_genes_20$DE_expression_series, collapse = '_')
  cis_20
  dim(de_genes_20)
  head(cis_20)
  library(dplyr)
  cis_20 <- cis_20 %>% left_join(de_genes_20,by=c("gene"="gene_symbol"))
  View(cis_20)
  
  
  
  cosmic_genes <- read.table(paste0(base_dir, 'biodatabase/oncogene_cosmic.txt'),sep='\t',header = T, check.names = F, stringsAsFactors = F)
  ref_cis_genes <- read.table(paste0(base_dir, 'biodatabase/cisplatin_resistance_genes.txt'),sep='\t',header = T, check.names = F, stringsAsFactors = F)
  # View(head(cosmic_genes))
  # View(head(ref_cis_genes))
  dim(cosmic_genes)
  dim(ref_cis_genes)
  cancer_genes_ls <- list()
  cosmic_genes_ls <- list()
  # rownames(pair_groups) <- pair_groups$desc
  minLogFC <- 0.25
  FDR_cutoff <- 0.01
  pValueThrs <- 0.05
  for(de in pair_groups$desc){
    signif_genes <- read.csv(paste0(pair_groups[de,'result_dir'],'signif_genes.csv'), check.names = F, stringsAsFactors=F)
    signif_genes <- signif_genes[abs(signif_genes$logFC)>minLogFC & signif_genes$FDR<FDR_cutoff  & signif_genes$PValue<pValueThrs,]
    signif_genes_cosmic <- signif_genes[signif_genes$gene_symbol %in% cosmic_genes$Gene_Symbol,]
    signif_genes_cis <- signif_genes[signif_genes$gene_symbol %in% ref_cis_genes$gene_symbol,]
    print(dim(signif_genes_cosmic))
    # fitness_genes <- summary(as.factor(signif_genes$classified_gene_dlp))
    nb_incis_cosmic <- sum(!is.na(signif_genes_cosmic$classified_gene_dlp))
    nb_intrans_cosmic <- sum(is.na(signif_genes_cosmic$classified_gene_dlp))
    # cancer_genes_ls[[paste0(datatag,'_',de)]] <- round(as.numeric(fitness_genes['TRUE'])*100/nrow(cancer_ref_genes_df),2)
    cosmic_genes_ls[[paste0(de,'_incis')]] <- round(nb_incis_cosmic*100/nrow(cosmic_genes),2)
    cosmic_genes_ls[[paste0(de,'_intrans')]] <- round(nb_intrans_cosmic*100/nrow(cosmic_genes),2)
    
    nb_incis_cis <- sum(!is.na(signif_genes_cis$classified_gene_dlp))
    nb_intrans_cis <- sum(is.na(signif_genes_cis$classified_gene_dlp))
    # cancer_genes_ls[[paste0(datatag,'_',de)]] <- round(as.numeric(fitness_genes['TRUE'])*100/nrow(cancer_ref_genes_df),2)
    cancer_genes_ls[[paste0(de,'_incis')]] <- round(nb_incis_cis*100/nrow(ref_cis_genes),2)
    cancer_genes_ls[[paste0(de,'_intrans')]] <- round(nb_intrans_cis*100/nrow(ref_cis_genes),2)
  }
  # summary(as.factor(pair_groups$title))
  # unique(pair_groups$title)
  plot_cisplatin_resistance_genes_pct(cancer_genes_ls)
  plot_cosmic_genes_pct(cosmic_genes_ls)
  
  
}

plot_cisplatin_resistance_genes_pct <- function(cancer_genes_ls){
  cis_df <- do.call(rbind, cancer_genes_ls)
  cis_df <- as.data.frame(cis_df)
  colnames(cis_df)[which(colnames(cis_df)=='V1')] <- 'cisplatin_resistance_gene_pct'
  dim(cis_df)
  View(head(cis_df))
  # colnames(cis_df) <- 'cosmic_genes_pct'
  
  # View(head(cis_df))
  level_pdx <- c('SA609_CISPLATIN','SA1035_CISPLATIN','SA535_CISPLATIN','SA535_CX5461')
  pdx <- lapply(strsplit(rownames(cis_df), "_"), function(x) {
    return(x[1])
  })
  cis_df$PDX <- as.character(pdx)
  unique(cis_df$PDX)
  cis_df$PDX <- ifelse(grepl('X',rownames(cis_df)),'SA535_CX5461',
                       ifelse(!grepl('X',rownames(cis_df)),paste0(cis_df$PDX,'_CISPLATIN'),cis_df$PDX))
  summary(as.factor(cis_df$PDX))
  # cis_df$PDX <- c(rep(level_pdx, c(12, 8, 14, 12)))
  cis_df$PDX <- factor(cis_df$PDX, levels=level_pdx)
  
  colnames(cancer_pct_df)[which(colnames(cancer_pct_df)=='V1')] <- 'Fitness_genes_pct'
  cis_df$desc <- as.character(rownames(cis_df)) #factor(rn, levels = rn)
  output_dir <- input_dir
  # View(head(cis_df))
  write.csv(cis_df, paste0(output_dir,'reference_cisplatin_resistance_gene_pct.csv'), quote=F, row.names = F)
  
  cis_df$gene_type <- ifelse(grepl("*incis",cis_df$desc),"incis",'intrans')
  # viz_pancancer_genes(cancer_pct_df, output_dir, tag="totaldata")
  summary(as.factor(cis_df$gene_type))
  
  cis_df$de_pair <- cis_df$desc
  # cancer_pct_df$de_pair <- gsub('_incis', '', grep("*incis",cancer_pct_df$de_pair, value=T))
  # cancer_pct_df$de_pair <- gsub('_intrans', '', grep("*intrans",cancer_pct_df$de_pair, value=T))
  cis_df$de_pair <- ifelse(grepl("*incis",cis_df$de_pair), gsub('_incis', '', cis_df$de_pair),
                           ifelse(grepl("*intrans",cis_df$de_pair), gsub('_intrans', '', cis_df$de_pair),cis_df$de_pair))
  
  cis_df$de_pair <- gsub('^(SA609_)|(SA1035_)|(SA535_)','',cis_df$de_pair)
  colorcode <- c('#FF3232','#40A0E0')
  ylabel <- ' (%) Ref cisplatin resistance genes'
  xlabel <- 'DE Analysis: Resistant versus Sensitive cells'
  plottitle <- 'Percentage DE genes in reference cisplatin resistance genes'
  View(head(cis_df))
  # cancer_pct_df$de_pair <- factor(cancer_pct_df$de_pair, levels = unique(cancer_pct_df$de_pair))
  cis_df$cisplatin_genes_pct
  # cancer_pct_df$series <- factor(cancer_pct_df$series, levels = unique(cancer_pct_df$series))
  plot_stack_barplot(cis_df, colorcode, xlabel, ylabel, plottitle, output_dir,'ref_cisplatin_resistance_genes_3series',
                     fa='gene_type', xa='de_pair', ya='cisplatin_resistance_gene_pct')
  
}  
plot_cosmic_genes_pct <- function(cosmic_genes_ls){
  cosmic_df <- do.call(rbind, cosmic_genes_ls)
  cosmic_df <- as.data.frame(cosmic_df)
  dim(cosmic_df)
  colnames(cosmic_df)[which(colnames(cosmic_df)=='V1')] <- 'cosmic_genes_pct'
  level_pdx <- c('SA609_CISPLATIN','SA1035_CISPLATIN','SA535_CISPLATIN','SA535_CX5461')
  pdx <- lapply(strsplit(rownames(cosmic_df), "_"), function(x) {
    return(x[1])
  })
  cosmic_df$PDX <- as.character(pdx)
  unique(cosmic_df$PDX)
  cosmic_df$PDX <- ifelse(grepl('X',rownames(cosmic_df)),'SA535_CX5461',
                          ifelse(!grepl('X',rownames(cosmic_df)),paste0(cosmic_df$PDX,'_CISPLATIN'),cosmic_df$PDX))
  summary(as.factor(cosmic_df$PDX))
  # cis_df$PDX <- c(rep(level_pdx, c(12, 8, 14, 12)))
  cosmic_df$PDX <- factor(cosmic_df$PDX, levels=level_pdx)
  
  
  cosmic_df$desc <- as.character(rownames(cosmic_df)) #factor(rn, levels = rn)
  output_dir <- input_dir
  View(head(cosmic_df))
  write.csv(cosmic_df, paste0(output_dir,'reference_cosmic_gene_pct.csv'), quote=F, row.names = F)
  cosmic_df$de_pair <- cosmic_df$desc
  # cancer_pct_df$de_pair <- gsub('_incis', '', grep("*incis",cancer_pct_df$de_pair, value=T))
  # cancer_pct_df$de_pair <- gsub('_intrans', '', grep("*intrans",cancer_pct_df$de_pair, value=T))
  cosmic_df$de_pair <- ifelse(grepl("*incis",cosmic_df$de_pair), gsub('_incis', '', cosmic_df$de_pair),
                              ifelse(grepl("*intrans",cosmic_df$de_pair), gsub('_intrans', '', cosmic_df$de_pair),cosmic_df$de_pair))
  cosmic_df$gene_type <- ifelse(grepl("incis",cosmic_df$desc),"incis",'intrans')
  cosmic_df$de_pair <- gsub('^(SA609_)|(SA1035_)|(SA535_)','',cosmic_df$de_pair)
  colorcode <- c('#FF3232','#40A0E0')
  ylabel <- '(%) Reference cosmic genes '
  xlabel <- 'DE Analysis: Resistant versus Sensitive cells'
  plottitle <- 'Percentage DE genes in COSMIC reference genes'
  
  # cancer_pct_df$de_pair <- factor(cancer_pct_df$de_pair, levels = unique(cancer_pct_df$de_pair))
  # cancer_pct_df$series <- factor(cancer_pct_df$series, levels = unique(cancer_pct_df$series))
  plot_stack_barplot(cosmic_df, colorcode, xlabel, ylabel, plottitle, output_dir,'ref_cosmic_genes_in_3series',
                     fa='gene_type', xa='de_pair', ya='cosmic_genes_pct')
  
}
plot_dlp_prevalence_v2 <- function(output_dir=NULL){
  output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_total_rna_v2/SA535-v6/'
  base_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/')
  cnv_535_fn <- paste0(base_dir, 'SA535_total_rna_v2/clonealign/whole_data/combined_clones/SA535_total_cnv_mat.csv')
  cnv_535 <- read.csv(cnv_535_fn, check.names = F, stringsAsFactors = F)
  head(cnv_535)
  cnv_1035_fn <- paste0(base_dir, 'SA1035_rna/clonealign/whole_data/SA1035_cnv_mat.csv')
  cnv_1035 <- read.csv(cnv_1035_fn, check.names = F, stringsAsFactors = F)
  
  cnv_609_fn <- paste0(base_dir, 'SA609_rna/added_segments/clonealign/whole_data/SA609_cnv_mat.csv')
  cnv_609 <- read.csv(cnv_609_fn, check.names = F, stringsAsFactors = F)
  
  
  de_df <- read.csv(paste0(base_dir, 'SA535_total_rna_v2/SA535-v6/comparisons_drug_res.csv'), check.names = F, stringsAsFactors = F)
  dim(de_df)  
  rownames(de_df) <- de_df$desc
  
  # de <- de_df$desc[8]
  cn_change <- list()
  for(de in de_df$desc){
    datatag <- de_df[de,'datatag']
    if(datatag=='SA609'){
      df_cnv <- cnv_609
    }else if(datatag=='SA535'){
      df_cnv <- cnv_535
    }else{
      df_cnv <- cnv_1035
    }
    df_cnv <- df_cnv %>% 
      dplyr::select(de_df[de,'clone1'],de_df[de,'clone2'])
    
    print(dim(df_cnv))
    df_cnv_var <- df_cnv[abs(df_cnv[,de_df[de,'clone1']]-df_cnv[,de_df[de,'clone2']])>=1,]
    print(df_cnv_var)
    # df_cnv_var <- df_cnv %>%
    #   dplyr::filter(apply(df_cnv, 1, var)>0)
    # dim(df_cnv)
    # dim(df_cnv_var)
    cn_change[[as.character(de)]] <- c(nrow(df_cnv_var),nrow(df_cnv))
  }
  
  cn_change_df <- as.data.frame(t(dplyr::bind_rows(cn_change)))
  colnames(cn_change_df) <- c('CN_change','total_CNs')
  cn_change_df$desc <- rownames(cn_change_df)
  write.csv(cn_change_df, paste0(output_dir,'cn_change.csv'), row.names = F, quote = F)
  View(cn_change_df)
  pdx <- lapply(strsplit(rownames(cn_change_df), "_"), function(x) {
    return(x[1])
  })
  cn_change_df$PDX <- as.character(pdx)
  cn_change_df$pct_change <- round(cn_change_df$CN_change/cn_change_df$total_CNs * 100,2)
  cn_change_df$PDX <- factor(cn_change_df$PDX, levels=c('SA609', 'SA1035','SA535'))
  p<-ggplot(data=cn_change_df, aes(x=desc, y=pct_change)) +
    geom_bar(stat="identity", fill="steelblue")+
    ggplot2::facet_grid(. ~ PDX, scales="free", space='free') +
    ggplot2::theme(legend.title = ggplot2::element_text(color="black", size=11),
                   legend.text = ggplot2::element_text(color="black", size=10),
                   plot.title = ggplot2::element_text(color="black", size=14, hjust = 0.5, face='bold'),
                   # legend.position = "none", 
                   # axis.line = element_blank(), 
                   panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), 
                   panel.border = ggplot2::element_blank(),
                   # axis.text.x =ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(color="black", size=10, hjust=0.5, angle = 90),
                   axis.text.y = ggplot2::element_text(color="black", size=8, hjust = 0.5),
                   axis.title.y = ggplot2::element_text(color="black", size=10, hjust = 0.5),
                   axis.title.x = ggplot2::element_text(color="black", size=12, hjust = 0.5, face='bold')
                   # axis.title.x = ggplot2::element_blank()
                   # axis.ticks = element_blank()
    )
  p
  
  
  
}  


plot_DE_genes_ggplot <- function(df, topGenes, capstr='', FDRcutoff=0.01, logFCcutoff=0.25, pValuecutoff=0.05,
                                 plttitle="A versus B", save_dir="",legendVisible=F,
                                 iscaption=TRUE, legend_verbose='none', save_plot=TRUE, 
                                 xl=NULL, yl=NULL){  
  # library(ggplot2)
  # library(ggrepel)
  #legend_verbose is none or 'right', 'left','bottom'
  # df <- de_genes
  # library(EnhancedVolcano)
  
  # colnames(df)[which(names(df) == "avg_logFC")] <- "log2FoldChange"
  # colnames(df)[which(names(df) == "p_val_adj")] <- "padj"
  # capstr <- paste0("FDR cutoff, ",FDRcutoff,"; logFC cutoff, ",logFCcutoff, "; nb genes signif, ",nbgenessig)
  # summary(as.factor(df$gene_type))
  Gene_Type=c('In_cis_Decrease_DownRegulated','In_cis_Decrease_UpRegulated',
              'In_cis_Increase_DownRegulated','In_cis_Increase_UpRegulated',
              'In_trans_DownRegulated','In_trans_UpRegulated'
  )
  gt <- data.frame(Gene_Type=Gene_Type,
                   gene_type=c('In cis Loss-Down','In cis Loss-Up',
                               'In cis Gain-Down','In cis Gain-Up',
                               'In trans Down','In trans Up'),
                   gene_type_classify=c('In cis positive tendency','In cis negative tendency',
                                        'In cis negative tendency','In cis positive tendency',
                                        'In trans','In trans'),
                   col_gt=c("#C3D7A4","#52854C",
                            "#FFDB6D","#D16103",
                            "#4E84C4","#293352"))
  
  df <- df %>% left_join(gt, by='Gene_Type')
  # unique(df$col_gt)
  # keyvals_colour <- as.character(df$col_gt)
  # names(keyvals_colour) <- df$gene_type
  # keyvals_colour <- gt$col_gt
  # names(keyvals_colour) <- gt$gene_type
  
  keyvals_colour <- c("In cis Loss-Down"="#C3D7A4","In cis Loss-Up"="#52854C",
                      "In cis Gain-Down"="#FFDB6D","In cis Gain-Up"="#D16103",
                      "In trans Down"="#4E84C4","In trans Up"="#293352")
  keyvals_colour_classify <- c("In cis positive tendency"="#6a0dad",
                               "In cis negative tendency"="#A8A8A8",
                               "In trans"="#000000")
  # keyvals_colour <- factor(keyvals_colour, levels = unique(keyvals_colour))
  # unique(keyvals_colour)
  # names(keyvals_colour[1:3])
  df <- df[abs(df$logFC)>logFCcutoff & df$FDR<FDRcutoff  & df$PValue<pValuecutoff,]
  maxLogFC <- 3.5
  # df$logFC <- ifelse(df$logFC>maxLogFC,maxLogFC,
  #                    ifelse(df$logFC<-maxLogFC,-maxLogFC,df$logFC))
  df$logFC <- sapply(df$logFC, function(x) replace(x, x > maxLogFC, maxLogFC))
  df$logFC <- sapply(df$logFC, function(x) replace(x, x < (-maxLogFC), -maxLogFC))
  
  st <- summary(df$gene_type)
  # keyvals_shape <- ifelse(df$is_fitness_gene==T, 3, 16)
  # names(keyvals_shape) <- ifelse(df$is_fitness_gene==T,'Fitness genes','Others')
  
  # df$gene_type <- factor(df$gene_type, levels = unique(df$gene_type))
  # if(capstr=='' & iscaption){
  #   # capstr <- paste0(capstr,'With abs(logFC)>0.25, FDR<0.01, pValue<0.05 \n')
  #   capstr <- paste0(capstr,names(st[1]),':',as.numeric(st[1]), ', ')
  #   capstr <- paste0(capstr,names(st[2]),':',as.numeric(st[2]), ' \n')
  #   capstr <- paste0(capstr,names(st[3]),':',as.numeric(st[3]), ', ')
  #   capstr <- paste0(capstr,names(st[4]),':',as.numeric(st[4]), ' \n')
  #   capstr <- paste0(capstr,names(st[5]),':',as.numeric(st[5]), ', ')
  #   capstr <- paste0(capstr,names(st[6]),':',as.numeric(st[6]), ' ')
  #   # for(i in rep(1:length(st),1)){
  #   #   # print(st[i])
  #   #   capstr <- paste0(capstr,names(st[i]),':',as.numeric(st[i]), ' ')
  #   # }
  # 
  # }
  df$mlog10FDR <- -log10(df$FDR)
  
  df$mlog10FDR <- sapply(df$mlog10FDR, function(x) replace(x, is.infinite(x), 300))
  # df$gt_alpha <- ifelse(df$gene_type=='In trans Up', 0.8,
  #                       ifelse(df$gene_type=='In cis Gain-Up',0.8,0.65))  
  # df$gt_alpha <- 0.8
  # df$gt_alpha <- ifelse(df$gene_type=='In trans Up' | df$gene_type=='In cis Gain-Up', 0.8, 0.9)  
  
  my_font <- "Helvetica"
  
  
  thesis_theme <- theme(  text=element_text(size = 8,family=my_font),
                          plot.title = element_text(color="black", size=11, hjust = 0, face = "bold", family=my_font),
                          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                          axis.text.y = element_text(size=7, hjust = 0.5, family=my_font),
                          axis.text.x = element_text(size=7, hjust = 0.5, family=my_font),
                          axis.title = element_text(size=8, hjust = 0.5, family=my_font),
                          plot.caption = element_text(size=8, hjust = 1, family=my_font),
                          legend.title = element_text(size=7, hjust = 0.5, family=my_font, angle=90),
                          legend.text = element_text(size=7, hjust = 0, family=my_font),
                          strip.text.x = element_text(color="black",size=9, family=my_font),
                          strip.text.y = element_text(color="black",size=9, family=my_font),
                          legend.spacing.x = unit(0.1, 'mm'),
                          legend.spacing.y = unit(0.1, 'mm'),
                          legend.key.height=unit(1,"line"),
                          # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          legend.position = legend_verbose)
  df_cis <- df %>%
    dplyr::filter(grepl('In cis',gene_type))
  
  topGenes <- get_top_genes(df_cis, minLogFC=0.25, nbtopup=10, nbtopdown=10)
  p_cis <- ggplot(df_cis, aes(x = logFC, y = mlog10FDR)) +
    geom_point(aes(color = gene_type), size=2.5, alpha=0.6) +
    scale_color_manual(values = keyvals_colour[1:4], name = "Gene Type") +
    thesis_theme #+ 
    # geom_text_repel(family = my_font,
    #                 data = df_cis[df_cis$gene_symbol %in% topGenes, ], 
    #                 aes(label = gene_symbol), size = 3.5, 
    #                 box.padding = unit(0.35, "lines"), 
    #                 point.padding = unit(0.3, "lines"),
    #                 max.overlaps = Inf,
    #                 min.segment.length = 0,  # draw segment lines, not matter how short they are)
    #                 color='black', segment.alpha = 0.2)
  p_cis <- p_cis + labs(x= bquote(~Log[2] ~ ' Fold Change'), 
                        y=bquote(~-Log[10]~italic(FDR)),title = paste0(plttitle,', ','in cis genes'))#caption = capstr
 
  
  if(!is.null(xl)){
    p_cis <- p_cis + xlim(xl[1],xl[2])
  }
  if(!is.null(yl)){
    p_cis <- p_cis + ylim(yl[1],yl[2])
  }
  
  
  df_trans <- df %>%
    dplyr::filter(grepl('In trans',gene_type))
  topGenes <- get_top_genes(df_trans, minLogFC=0.25, nbtopup=10, nbtopdown=10)
  p_trans <- ggplot(df_trans, aes(x = logFC, y = mlog10FDR)) +
    geom_point(aes(color = gene_type), size=2.5, alpha=0.6) +
    scale_color_manual(values = keyvals_colour[5:6], name = "Gene Type") +
    thesis_theme #+ 
    # geom_text_repel(family = my_font,
    #                 data = df_trans[df_trans$gene_symbol %in% topGenes, ], 
    #                 aes(label = gene_symbol), size = 3.5, 
    #                 box.padding = unit(0.35, "lines"), 
    #                 point.padding = unit(0.3, "lines"),
    #                 max.overlaps = Inf,
    #                 min.segment.length = 0,  # draw segment lines, not matter how short they are)
    #                 color='black', segment.alpha = 0.2)
  p_trans <- p_trans + labs(x= bquote(~Log[2] ~ ' Fold Change'), 
                            y=bquote(~-Log[10]~italic(FDR)),title = paste0(plttitle,', ','in trans genes'))#caption = capstr
  
  if(!is.null(xl)){
    p_trans <- p_trans + xlim(xl[1],xl[2])
  }
  if(!is.null(yl)){
    p_trans <- p_trans + ylim(yl[1],yl[2])
  }
  
  # if(legend_verbose!='none'){
  #   return(p)
  # }
  
  # xdens <- cowplot::axis_canvas(p, axis = "x") +
  #   stat_density(data = df, aes(x = logFC, group = gene_type_classify,color =gene_type_classify),
  #                alpha = 1, size = 1, position="identity",geom="line") +
  #   scale_color_manual(values = keyvals_colour_classify, name = "Gene Classify") 
  
  
  # xdens <- cowplot::axis_canvas(p, axis = "x") +
  #          # geom_jitter(data = df, aes(x = logFC, fill = gene_type_classify), 
  #          #             shape=16, position=position_jitter(0.2)) + 
  #          geom_boxplot(data = df, aes(x = logFC, color = gene_type_classify),
  #                alpha = 0.6, size = 0.3) +
  #         # scale_color_grey()
  #         scale_color_manual(values = keyvals_colour_classify, name = "Gene Classify")
  #         
  # 
  # p_main <- cowplot::insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
  # 
  
  # Just to get legend
  # t <- ggplot(df, aes(x = logFC, color = gene_type_classify)) +
  #      geom_boxplot(alpha = 0.6, size = 0.3) +
  #      scale_color_manual(values = keyvals_colour_classify, name = "Gene Classify") + 
  #      theme(text=element_text(size = 8,family=my_font),
  #           legend.title = element_text(size=7, hjust = 0.5, family=my_font),
  #           legend.text = element_text(size=7, hjust = 0, family=my_font),
  #           strip.text.x = element_text(color="black",size=9, family=my_font),
  #           strip.text.y = element_text(color="black",size=9, family=my_font),
  #           legend.spacing.x = unit(0.1, 'mm'),
  #           legend.spacing.y = unit(0.1, 'mm'),
  #           legend.key.height=unit(1,"line"))
  # lg <- cowplot::get_legend(t)
  # plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  
  # blank_theme <- theme_minimal()+
  #   theme(
  #     axis.title.x = element_blank(),
  #     axis.title.y = element_blank(),
  #     panel.border = element_blank(),
  #     panel.grid=element_blank(),
  #     axis.ticks = element_blank(),
  #     plot.title=element_text(size=14, face="bold")
  #   )
  # prop_df <- df %>%
  #   dplyr::group_by(gene_type)%>%
  #   dplyr::summarise(count=n())%>%
  #   as.data.frame()%>%
  #   dplyr::rename(label=gene_type)
  
  # For plotting
  shot_lbs <- c('LD','LU','GD','GU','Down','Up')
  names(shot_lbs) <- names(keyvals_colour)
  prop_cis <- df %>%
    dplyr::filter(grepl('In cis',gene_type)) %>%
    dplyr::group_by(gene_type)%>%
    dplyr::summarise(count=n())%>%
    as.data.frame()#%>%
  # dplyr::rename(label1=gene_type)
  prop_cis$label <- shot_lbs[prop_cis$gene_type]
  
  prop_trans <- df %>%
    dplyr::filter(grepl('In trans',gene_type)) %>%
    dplyr::group_by(gene_type)%>%
    dplyr::summarise(count=n())%>%
    as.data.frame()#%>%
  # dplyr::rename(label1=gene_type)
  prop_trans$label <- shot_lbs[prop_trans$gene_type]
  # dim(prop_cis)
  names(keyvals_colour) <- shot_lbs[names(keyvals_colour)]
  pcis <- viz_proportion(prop_cis, plot_label=F, keyvals_colour)
  ptrans <- viz_proportion(prop_trans, plot_label=F, keyvals_colour)
  # p1
  
  # p1 <- ggplot(df, aes(x="", y=gene_type, fill=gene_type))+
  #   geom_bar(width = 1, stat = "identity", alpha=0.8) +
  #   coord_polar("y", start=0) + 
  #   scale_fill_manual(values=keyvals_colour) + 
  #   blank_theme + 
  #   theme(axis.text.x=element_blank(), legend.position = "none")
  # 
  # caption = capstr
  # plg <- ggplot() +
  #   annotate("text", x = 8, y = 20, color="black", size=4.5, family=my_font, label = capstr) +
  #   theme_void()
  # p2
  # p1 <- cowplot::plot_grid(p2, p1, ncol=2, rel_widths = c(2,1))
  # p_total <- cowplot::plot_grid(p, p1, rel_heights = c(1,0.2), ncol=1)
  # p_total
  # tag <- ifelse(legend_verbose=='none','','_with_legend')
  # png(paste0(save_dir,"DE_",gsub(' ','_',plttitle),tag,".png"), height = 2*650, width=2*500, res = 2*72)
  # print(p_total)
  # dev.off()
  
  
  # if(save_plot){
  #   plttitle <- gsub(':','_',plttitle)
  #   plttitle <- gsub(' ','_',plttitle)
  #   saveRDS(p_total, file=paste0(save_dir,"DE_",plttitle,".rds"))
  #   ggsave(paste0(save_dir,"DE_",gsub(' ','_',plttitle),tag,".pdf"),
  #          plot = p_total,
  #          height = 5.5,
  #          width = 6.5,
  #          useDingbats=F)
  # }
  return(list(p_vol_cis=p_cis, p_vol_trans=p_trans, p_prop_cis=pcis, p_prop_trans=ptrans))
}

viz_legend <- function(df, keyvals_colour){
  p1 <- ggplot(df, aes(x = logFC, y = mlog10FDR)) +
    geom_point(aes(color = gene_type), size=2.5,alpha=0.8) +
    scale_color_manual(values = keyvals_colour, name = "Gene Type") + 
    theme(legend.position ='bottom')
  p1 <- p1 + guides(color = guide_legend(override.aes = list(size=4, nrow=3)))
  lg <- cowplot::get_legend(p1)
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  return(plg)
}
# Donut plot
# Need count, label
viz_proportion <- function(data, plot_label=F, cols=NULL){
  data <- data %>%
    dplyr::filter(count>0)
  # if(is.null(cols)){
  #   cols blah blah
  # }
  
  # Compute percentages
  # data$label <- factor(data$label, levels = names(cols))
  rownames(data) <- data$label
  data <- data[names(cols)[names(cols) %in% as.character(data$label)],]
  data$fraction <- data$count / sum(data$count)
  
  # Compute the cumulative percentages (top of each rectangle)
  data$ymax <- cumsum(data$fraction)
  
  # Compute the bottom of each rectangle
  data$ymin <- c(0, head(data$ymax, n=-1))
  
  # Compute label position
  data$labelPosition <- (data$ymax + data$ymin) / 2
  my_font <- "Helvetica"
  # Compute a good label
  # data$label <- data$category
  # data$label <- c('1Rx','2Rx','3Rx','1RxH')
  cols_use <- cols[as.character(unique(data$label))]
  data$label_desc <- paste0(data$label, "\n (", data$count,')') #\n: N=
  if(plot_label==F){
    p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=label)) +
      geom_rect() +
      # geom_label( x=3.5, aes(y=labelPosition, label=label), size=3.1, label.size=0) +
      scale_fill_manual(values = cols_use)   
  }else{
    p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=label)) +
      geom_rect() +
      # geom_label(x=5.5, aes(y=labelPosition, label=label), size=3, label.size=0) +
      geom_text( x=5, aes(y=labelPosition, label=label_desc, color=label), size=3.5, family=my_font) +
      scale_fill_manual(values = cols_use) + 
      scale_color_manual(values = cols_use) 
  }
  p <- p + coord_polar(theta="y") +
    xlim(c(1.5, 5.5)) +
    theme_void() +
    theme(legend.position = "none",
          text=element_text(family=my_font))
  # plot.margin = margin(0, 0, 0, 0, "cm"),
  # axis.ticks.length = unit(0, "null")
  # panel.background = element_rect(fill='transparent'),
  # plot.background = element_rect(fill='transparent', color=NA),
  # panel.grid.major = element_blank(),
  # panel.grid.minor = element_blank()
  # p
  return(p)
  
}


viz_pancancer_genes <- function(cancer_pct_df, output_dir, tag=""){
  p <- ggplot(data=cancer_pct_df, aes(x=PDX, y=Fitness_genes_pct, fill=series)) +
    geom_bar(stat="identity", width=0.7)+ #, fill="white"
    geom_text(aes(label=Fitness_genes_pct), vjust=-0.3, size=3.5, color="black")+
    theme_minimal()
  p <- p + scale_fill_manual(values=c("#30F50A", "#8DD87F", "#1A6A0B", "#560C56"))
  
  p <- p + labs(x='PDX drug treatment', y="% Fitness_genes", title='Percentage genes in ADAM PanCancer Core-Fitness genes')
  p <- p + theme(legend.title = element_text(size=10), 
                 legend.text = element_text(size=9),
                 plot.title = element_text(color="black", size=13, hjust = 0.5),
                 # legend.position = "none", 
                 # axis.line = element_blank(), 
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 # panel.border = element_blank(),
                 axis.text.x = element_text(size=9, angle = 90, color="black"),
                 #axis.title = element_blank(),
                 axis.ticks.x = element_blank()
  )
  # p
  png(paste0(output_dir,tag,"_fitness_genes_stat.png"), height = 2*400, width=2*620, res = 2*72)
  print(p)
  dev.off()
}

plot_reference_set_proportion_total_v2 <- function(stat, save_dir, verbose=F){
  
  prop_plt_ls <- list()
  prop_plt <- plot_reference_set_proportion_v2(stat_tmp, save_dir, verbose=verbose)
  
  return(prop_plt)
  
}

plot_reference_set_proportion_total <- function(stat, save_dir, verbose=F){
  datatags <- c('SA609','SA535','SA1035')
  prop_plt_ls <- list()
  for(tag in datatags){
    stat_tmp <- stat %>%
      dplyr::filter(datatag==tag & comp_type=='treated_vs_untreated')
    # dim(stat_tmp)
    # prop_plt <- plot_reference_set_proportion(stat_tmp, save_dir, verbose=verbose)
    prop_plt <- plot_reference_set_proportion_v2(stat_tmp, save_dir, verbose=verbose)
    prop_plt_ls[[tag]] <- prop_plt$prevalence_plt
    # prop_plt$prevalence_plt
  }
  prop_plt_ls$plg <- prop_plt$plg
  return(prop_plt_ls)
  
}

plot_pathway_custom_sets_total <- function(pw_stat, pair_groups, save_dir, verbose=F){
  datatags <- c('SA609','SA535','SA1035')
  pathway_plt_ls <- list()
  for(tag in datatags){
    groups_use <- pair_groups %>%
      dplyr::filter(datatag==tag & comp_type=='treated_vs_untreated')
    # groups_use$order
    dim(groups_use)
    stat_tmp <- pw_stat %>%
      dplyr::filter(file_header %in% groups_use$file_header)
    
    # groups_use <- groups_use %>% left_join(stat_tmp, by='file_header')
    stat_tmp <- stat_tmp %>% left_join(groups_use, by='file_header')
    # dim(groups_use)
    # groups_use$padj
    # View(groups_use)
    # prop_plt <- viz_reference_pathways(groups_use, verbose=verbose)
    prop_plt <- viz_reference_pathways(stat_tmp, verbose=verbose)
    pathway_plt_ls[[tag]] <- prop_plt$pathway_sig_plt
    # prop_plt$prevalence_plt
  }
  pathway_plt_ls$plg <- prop_plt$plg
  return(pathway_plt_ls)
  
}

plot_reference_set_proportion_v2 <- function(stat, save_dir, verbose=FALSE, plttype=NULL){
  alias <- data.frame(ref_set=c('CF','BroadSanger','cosmic','CR'), 
                      ref_gene_set=c('CF','BS Essential','COSMIC','CR')) #'Core-Fitness(CF)', 'Cispl-Res(CR)'
  
  stat <- stat %>% 
    dplyr::filter(ref_set %in% c('CF','CR'))%>% 
    left_join(alias, by=c('ref_set'))
  # unique(stat$ref_gene_set)
  # sum(is.na(stat$ref_gene_set))
  my_font <- "Helvetica"
  # colorcode <- c('#FF3232','#40A0E0')
  colorcode <- c("#D2691E","blue2")
  names(colorcode) <- c("In cis","In trans")
  stat$gene_type <- ifelse(stat$gene_type=='In_cis','In cis','In trans')
  stat$proportion_size <- ifelse(stat$pct_genes==0, 0.01, log2(stat$pct_genes)) #0.01 is just for illustration
  stat$pct_genes <- ifelse(stat$pct_genes==0, 0.01, stat$pct_genes)
  stat <- stat %>%
    dplyr::arrange(-order)
  
  if(verbose==F){
    axisY <- ggplot2::element_blank()
  }else{
    axisY <- ggplot2::element_text(size=9, hjust = 0.5, family=my_font)
  }
  
  # stat$labels_detailed <- factor(stat$labels_detailed, levels = unique(stat$labels_detailed))
  p <- ggplot(stat, aes(x = gene_type, y = order)) + 
    geom_point(aes(color = gene_type, size=proportion_size), alpha=1) + #, size=2*log10(pct_genes)size=4.5
    scale_color_manual(values = colorcode) +
    scale_y_discrete(drop=T)
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    # scale_x_continuous(trans='log10') + 
    # facet_grid(series ~ ref_gene_set, scales="free_y", space='free',drop=T)  # PDX ~ ref_gene, . ~ PDX
    # geom_hline(yintercept = 1:6, col = "#e2e2e2") +
    # geom_text(aes(label=celltype_desc))+
    # annotate('text', x = df$success, y = df$index, label = df$success, size=3, colour = col[df$gender])+
  if(!is.null(plttype)){
    p <- p + facet_grid(series ~ ref_gene_set, scales="free_y", space='free',drop=T)
  }else{
    p <- p + facet_grid(. ~ ref_gene_set, scales="free_y", space='free',drop=T)
  }
    # scale_color_manual(values = col) +
  p <- p + 
    ct_theme
    # theme_bw() + 
    # theme(strip.text.x = element_text(size=11, color="black", family=my_font),
    #       strip.text.y = element_blank(),
    #       strip.background = element_blank(),
    #       legend.position = "none",
    #       # panel.grid.major.x = element_blank(),
    #       # panel.grid.minor.x = element_blank(),
    #       # axis.ticks.y = element_blank(),
    #       axis.text.y  = axisY,
    #       text = element_text(size = 10, hjust = 0.5, family=my_font),
    #       axis.text.x = element_text(size=11, hjust = 0.5, family=my_font),  #, angle = 90
    #       # axis.text.y = element_text(size=9, hjust = 0.5, family=my_font),
    #       plot.title = element_text(size=11, face="bold", hjust=0, family=my_font),
    #       axis.title.y = element_blank(),
    #       axis.title.x = element_text(size=10, hjust = 0.5, family=my_font))
  p <- p + labs(x=NULL,y=NULL, title = NULL) #'Custom gene set proportion (%)'
  # p
  # geom_text_repel(data = stat, aes(label = pct_gene), nudge_x=0, nudge_y=0, max.overlaps = Inf,
  #                  min.segment.length = 0, segment.size=0.2,
  #                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines") )
  # p
  # p <- p + guides(color=guide_legend(title="Gene Type ", override.aes = list(shape = 0, size=6))) #size = FALSE,
  
  
  
  p1 <- ggplot(stat, aes(x = pct_genes, y = labels_detailed)) + 
    geom_point(aes(color = gene_type), size=6.5, alpha=1) + #, size=2*log10(pct_genes)size=4.5
    scale_color_manual(name = 'Gene Type', values = colorcode, labels = c("In Cis", "In Trans")) + 
    theme(legend.position ='bottom')
  p1 <- p1 + guides(fill = guide_legend(override.aes = list(shape = 0, size=8, nrow = 1)))
  lg <- cowplot::get_legend(p1)
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  # plg
  res <- list(prevalence_plt=p, plg=plg)
  
  # saveRDS(res, paste0(save_dir,'proportion_ref_genes.rds'))
  # ptotal <- cowplot::plot_grid(p, plg, ncol = 1, rel_heights = c(15,1.3))
  # png(paste0(save_dir,"Fig2_part21.png"), height = 2*600, width=2*520, res = 2*72)
  # print(ptotal)
  # dev.off()
  
  return(res)
}

plot_reference_set_proportion <- function(stat, save_dir, verbose=FALSE){
  
  # pair_groups: for order of comparisons
  # alias <- data.frame(ref_set=c('CoreFitness','BroadSanger','cosmic','cisplatin_resistance'), 
  #                     ref_gene=c('Core Fitness (CF)','BS Essential','COSMIC','Cispl Resistance (CR)'))
  alias <- data.frame(ref_set=c('CF','BroadSanger','cosmic','CR'), 
                      ref_gene_set=c('Core Fitness (CF)','BS Essential','COSMIC','Cispl Resistance (CR)'))
  # stat <- stat %>% pivot_longer(cols = starts_with("in"),
  #                               names_to = 'Gene_Type', values_to = 'pct_gene')
  # dim(stat)
  
  stat <- stat %>% 
    dplyr::filter(ref_set %in% c('CF','CR'))%>% 
    left_join(alias, by=c('ref_set'))
  # unique(stat$ref_gene_set)
  # sum(is.na(stat$ref_gene_set))
  my_font <- "Helvetica"
  colorcode <- c('#FF3232','#40A0E0')
  names(colorcode) <- c('In_cis','In_trans')
  
  stat$proportion_size <- ifelse(stat$pct_genes==0, 0.01, log10(stat$pct_genes)) #0.01 is just for illustration
  stat$pct_genes <- ifelse(stat$pct_genes==0, 0.01, stat$pct_genes)
  stat <- stat %>%
    dplyr::arrange(-order)
  
  if(verbose==F){
    axisY <- ggplot2::element_blank()
  }else{
    axisY <- ggplot2::element_text(size=9, hjust = 0.5, family=my_font)
  }
  stat$labels_detailed <- factor(stat$labels_detailed, levels = unique(stat$labels_detailed))
  p <- ggplot(stat, aes(x = pct_genes, y = labels_detailed)) + 
    geom_point(aes(color = gene_type, size=proportion_size), alpha=0.9) + #, size=2*log10(pct_genes)size=4.5
    scale_color_manual(values = colorcode) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    # scale_x_continuous(trans='log10') + 
    facet_grid(~ ref_gene_set, scales="free_y", space='free',drop=T) + # PDX ~ ref_gene, . ~ PDX
    # geom_hline(yintercept = 1:6, col = "#e2e2e2") +
    # geom_text(aes(label=celltype_desc))+
    # annotate('text', x = df$success, y = df$index, label = df$success, size=3, colour = col[df$gender])+
    
    # scale_color_manual(values = col) +
    theme_bw() + 
    theme(strip.text = element_text(size=9, color="black", family=my_font),
          strip.background = element_blank(),
          legend.position = "none",
          # panel.grid.major.x = element_blank(),
          # panel.grid.minor.x = element_blank(),
          # axis.ticks.y = element_blank(),
          axis.text.y  = axisY,
          text = element_text(size = 7, hjust = 0.5, family=my_font),
          axis.text.x = element_text(size=9, hjust = 0.5, family=my_font),  #, angle = 90
          # axis.text.y = element_text(size=9, hjust = 0.5, family=my_font),
          plot.title = element_text(size=10, face="bold", hjust=0, family=my_font),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size=10, hjust = 0.5, family=my_font))
  p <- p + labs(x='Custom gene set proportion (%)',y=NULL, title = NULL)
  # p
  # geom_text_repel(data = stat, aes(label = pct_gene), nudge_x=0, nudge_y=0, max.overlaps = Inf,
  #                  min.segment.length = 0, segment.size=0.2,
  #                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines") )
  # p
  # p <- p + guides(color=guide_legend(title="Gene Type ", override.aes = list(shape = 0, size=6))) #size = FALSE,
  
  
  
  p1 <- ggplot(stat, aes(x = pct_genes, y = labels_detailed)) + 
    geom_point(aes(color = gene_type), size=6.5, alpha=1) + #, size=2*log10(pct_genes)size=4.5
    scale_color_manual(name = 'Gene Type', values = colorcode, labels = c("In Cis", "In Trans")) + 
    theme(legend.position ='bottom')
  p1 <- p1 + guides(fill = guide_legend(override.aes = list(shape = 0, size=8, nrow = 1)))
  lg <- cowplot::get_legend(p1)
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  # plg
  res <- list(prevalence_plt=p, plg=plg)
  
  # saveRDS(res, paste0(save_dir,'proportion_ref_genes.rds'))
  # ptotal <- cowplot::plot_grid(p, plg, ncol = 1, rel_heights = c(15,1.3))
  # png(paste0(save_dir,"Fig2_part21.png"), height = 2*600, width=2*520, res = 2*72)
  # print(ptotal)
  # dev.off()
  
  return(res)
}
# pathway_stat <- groups_use
viz_reference_pathways <- function(pathway_stat, verbose=F, plttype=NULL){
  # library(viridis)
  my_font <- "Helvetica"
  if(verbose==F){
    axisY <- ggplot2::element_blank()
  }else{
    axisY <- ggplot2::element_text(size=9, hjust = 0.5, family=my_font)
  }
  # pathway_stat$ref_gene <- ifelse(pathway_stat$ref_gene=='Cisplt Resistance','Cisplt Res',
  #                                 pathway_stat$ref_gene)
  # pathway_stat <- pathway_stat %>%
  #   dplyr::mutate(padj=replace(padj, padj>=0.05, 0))
  
  # pathway_stat$cell_size <- ifelse(pathway_stat$padj>0, log10(pathway_stat$size), NA)
  # pdx_level = c('SA501','SA530','SA604','SA609','SA535','SA1035')

  pathway_stat <- pathway_stat %>% 
    dplyr::mutate(pathway = case_when(pathway=='CISPLATIN_RESISTANCE' ~ 'CR',
                                      pathway=='CORE_FITNESS' ~ 'CF'),
                  cell_size = case_when(padj<0.05 ~ 2L,
                                        padj>=0.05 ~ NA_integer_))%>% 
    dplyr::arrange(-order)
  pathway_stat$pathway <- factor(pathway_stat$pathway, levels = c('CF','CR'))
  # pathway_stat$order
  
  # pathway_stat$type <- 'Stat Level'
  # pathway_stat$labels_detailed <- factor(pathway_stat$labels_detailed, levels=unique(pathway_stat$labels_detailed))
  p <- ggplot(pathway_stat, aes(x=pathway, y=order, color=padj)) +
    geom_point(aes(size=cell_size), shape=8) +
    viridis::scale_color_viridis(discrete=F, alpha=0.8, limits = c(min(pathway_stat$padj),0.05))+
    scale_y_discrete(drop=T)
    # facet_grid(PDX ~ type, scales="free_y", space='free',drop=T) + # PDX ~ ref_gene, . ~ PDX
    # theme_bw() + 
  if(!is.null(plttype)){
    p <- p + facet_grid(series ~ ., scales="free_y", space='free', drop=T)
  }  
  p <- p + 
    ct_theme
       # theme(strip.text = element_text(size=11, color="black", family=my_font),
       #    # strip.text = element_blank(),
       #    strip.background = element_blank(),
       #    legend.position = "none",
       #    # panel.grid.major.x = element_blank(),
       #    # panel.grid.minor.x = element_blank(),
       #    # axis.ticks.y = element_blank(),
       #    axis.text.y = axisY,
       #    # axis.text.y = element_text(size=9, hjust = 0.5, family=my_font, angle = 90),
       #    text = element_text(size = 11, hjust = 0.5, family=my_font),
       #    axis.text.x = element_blank(),
       #    # axis.text.x = element_text(size=11, hjust = 0.5, family=my_font, angle = 90, color="black"),  #, angle = 90
       #    # axis.text.y = element_text(size=9, hjust = 0.5, family=my_font),
       #    plot.title = element_text(size=10, hjust=0.5, family=my_font),
       #    # axis.title.x = element_text(size=9, hjust = 0.5, family=my_font),
       #    axis.title.x = element_blank(),
       #    axis.title.y = element_blank(),
       #    # plot.background = element_rect(fill = "white", colour = "white"),
       #    panel.background = element_rect(fill = 'white', colour = 'white')) 
  p <- p + labs(title=NULL, x=NULL, y=NULL) #'Custom Set''Significant'
  # p <- p + guides(color=guide_legend(title="Sigf padj")) #size = FALSE,
  # p
  p1 <- ggplot(pathway_stat, aes(x=pathway, y=order, color=padj)) +
    geom_point(size=5) +
    viridis::scale_color_viridis(discrete=F, alpha=0.8, name = 'P-adjusted  ', limits = c(min(pathway_stat$padj),0.05)) +
    theme(legend.position ='right',
          legend.key.width = unit(0.3, 'cm'),
          legend.key.height = unit(0.8, 'cm'), 
          legend.text = element_text(size = 8, hjust = 0.5, family=my_font))
  # p1 <- p1 + guides(fill = guide_legend(override.aes = list(shape = 0, size=8, nrow = 1)))
  
  lg <- cowplot::get_legend(p1)
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  # plg
  res_pathway <- list(pathway_sig_plt=p, plg=plg)
  # saveRDS(res_pathway, paste0(save_dir,'pathway_ref_plt.rds'))
  # ptotal <- cowplot::plot_grid(p, plg, ncol = 1, rel_heights = c(15,1))
  # png(paste0(save_dir,"Fig2_part22.png"), height = 2*600, width=2*150, res = 2*72)
  # print(ptotal)
  # dev.off()
  
  return(res_pathway)
  
}
