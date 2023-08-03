
# SA609
# script_dir <- '/home/htran/Projects/farhia_project/rnaseq/trajectory_analysis/'
script_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/scripts/'
source(paste0(script_dir, "trajectory_analysis/slingshot_utils.R"))
source(paste0(script_dir, "trajectory_analysis/tradeseq_utils.R"))
library(grid)
# output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/'
# results_dir <- paste0(output_dir,'figs_v3/')
# my_font <- "Helvetica"
# datatag <- 'SA609'
output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/'
results_dir <- paste0(output_dir,'figs_v5/')
save_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/figs_v5/"
my_font <- "Helvetica"
datatag <- 'SA609'

plot_trajectory <- function(){
  ## Trajectory plot, old version
  ## See output of res <- plot_all_lingeages(sce, crv_umap_embed, output_dir, datatag)
  # l1 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage1_SA609.rds'))
  # l1 <- l1 + theme(plot.title = element_text(color="black", size=12, hjust = 0.5, family=my_font))
  # 
  # l2 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage2_SA609.rds'))
  # l2 <- l2 + theme(plot.title = element_text(color="black", size=12, hjust = 0.5, family=my_font))
  # 
  # l3 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage3_SA609.rds'))
  # l3 <- l3 + theme(plot.title = element_text(color="black", size=12, hjust = 0.5, family=my_font))
  # 
  # total <- readRDS(paste0(results_dir,'ts_slingshot_out_wholedataset_SA609.rds'))
  # total <- total + theme(plot.title = element_text(color="black", size=14, hjust = 0, family=my_font, face="bold"))
  # ts_color <- readRDS(paste0(results_dir,'treatment_plt.rds'))
  # # ts_color <- lg2
  # 
  # 
  # p1 <- cowplot::plot_grid(l3, l2, l1, ncol=1)
  # p2 <- cowplot::plot_grid(total, ts_color, ncol=1, rel_heights = c(4,1))+
  #   theme(plot.background = element_rect(fill = "white", colour = "white"))
  # ptotal <- cowplot::plot_grid(p2, p1, rel_widths = c(2,1), nrow = 1)
  # ptotal
  # dev.off()
  # png(paste0(results_dir,"slingshot_out_summary_",datatag,"_part11.png"), height = 2*500, width=2*700,res = 2*72)
  # print(ptotal)
  # dev.off()
  
  ## Trajectory plot
  ## See output of res <- plot_all_lingeages(sce, crv_umap_embed, output_dir, datatag)
  l1 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage1_SA609.rds'))
  
  # p <- readRDS(paste0(results_dir,"ts_slingshot_out_",lid,"_",datatag,".rds"))
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
  total <- total + theme(plot.title = element_text(color="black", size=16, hjust = 0, family=my_font, face="bold"))
                         # axis.title = element_text(color="black", size=10, hjust = 0, family=my_font))
  ts_color <- readRDS(paste0(results_dir,'treatment_plt.rds'))
  # ts_color <- lg2
  # dev.off()
  ts_color <- readRDS(paste0(results_dir, '_plt.rds'))
  
  p1 <- cowplot::plot_grid(l3, l2, l1, ncol=1)
  p_plg <- cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_plot(res$hm_plg_exp), 
                              cowplot::ggdraw() + cowplot::draw_plot(res$hm_plg_cistrans), 
                              cowplot::ggdraw() + cowplot::draw_plot(res$hm_plg_chromatin),
                              ncol=3, rel_widths = c(1, 1, 1),
                              hjust = -0.1,label_fontface="plain") #+ 
    # theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  # p2 <- cowplot::plot_grid(total, ts_color, ncol=1, rel_heights = c(4,1))+
  #   theme(plot.background = element_rect(fill = "white", colour = "white"))
  p2 <- cowplot::plot_grid(total, p_plg, ncol=1, rel_heights = c(3.6,1))+
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  ptotal_traj_609 <- cowplot::plot_grid(p2, p1, rel_widths = c(2,1), nrow = 1)
  ptotal_traj_609
  
  
  
  ## SA609 most updated version
  p <- cowplot::plot_grid(l3, l2, l1, 
                          NULL, nrow=2) + 
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  p_total <- cowplot::plot_grid(ts_color, p, rel_heights = c(1, 8), ncol=1)+
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  # ggsave(paste0(save_dir,"hm_gene_modules_",datatag, ".png"),
  #        plot = ptotal_traj_609,
  #        height = 4,
  #        width = 5,
  #        # useDingbats=F,
  #        dpi=150)
  save_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/figs/'
  ggsave(paste0(save_dir,"umaps_lineages.png"),
         plot = p,
         height = 4,
         width = 5,
         # useDingbats=F,
         dpi=150)
  
  save_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/figs/'
  ggsave(paste0(save_dir,"umaps_lineages_SA609.png"),
         plot = p,
         height = 8,
         width = 10,
         # useDingbats=F,
         dpi=250)
  
  save_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/figs/'
  ggsave(paste0(save_dir,"umaps_lineages_SA535.png"),
         plot = p,
         height = 9,
         width = 12,
         # useDingbats=F,
         dpi=250)
  
  save_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/figs/'
  ggsave(paste0(save_dir,"umaps_lineages_SA535_plg.svg"),
         plot = lg2,
         height = 1.5,
         width = 7,
         # useDingbats=F,
         dpi=250)
  
  # ggsave(paste0(save_dir,"hm_gene_modules_",datatag, "_lg.png"),
  #        plot = ts_color,
  #        height = 3.5,
  #        width = 3.5,
  #        # useDingbats=F,
  #        dpi=150)
  
  return(ptotal_traj_609)
}

ptotal_traj_609 <- plot_trajectory()

plot_heatmap_gene_modules <- function(){
  # output_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/'
  # obs_genes_df1 <- data.table::fread(paste0(output_dir,"SA609_total_genes_modules_act_repr_trans_08_Dec.csv.gz")) %>% as.data.frame()
  # save_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/"
  # obs_genes_df <- data.table::fread(paste0(save_dir,"obs_genes_hm.csv.gz")) %>% as.data.frame()
  # obs_genes_df <- data.table::fread(paste0(output_dir,"SA609/obs_genes_hm.csv.gz")) %>% as.data.frame()
  # obs_genes_df$gene_type
  # dim(obs_genes_df)
  # summary(as.factor(obs_genes_df$gene_type))
  # dim(obs_genes_df1)
  # summary(as.factor(obs_genes_df1$gene_type_module))
  # ## Heatmap average expression plot
  ## Extracting heatmap of average genes expression x lineages 
  save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/'
  # genes_df <- data.table::fread(paste0(save_dir, 'SA609_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
  # dim(genes_df)
  # colnames(genes_df)
  # genes_df <- genes_df %>%
  #   dplyr::select(-gene_type) %>%
  #   dplyr::rename(gene_type=gene_type_module)
  # summary(as.factor(genes_df$gene_type))
  datatag <- 'SA609'
  plttitle <- 'Activated Repressed Transient genes'
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/')
  # ts_sce <- readRDS(paste0(save_dir,'tradeSeq_3000/', "fitGAM_out.rds"))
  save_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/"
  ## phm <- viz_heatmap(ts_sce, genes_df, save_dir, datatag, plttitle)
  # output_dir <- save_dir # old version
  output_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/SA609/'
  exp_mtx <- data.table::fread(paste0(output_dir,"mtx_hm.csv.gz")) %>% as.data.frame()
  obs_genes_df <- data.table::fread(paste0(output_dir,"obs_genes_hm.csv.gz")) %>% as.data.frame()
  obs_cells_df <- data.table::fread(paste0(output_dir,"obs_cells_hm.csv.gz")) %>% as.data.frame()
  rownames(exp_mtx) <- exp_mtx$ens_gene_id
  exp_mtx$ens_gene_id <- NULL
  input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/cis_trans_lineages/')
  # input_dir <- 'where to keep cistrans genes' // TO DO
  cistrans_anno <- data.table::fread(paste0(input_dir, datatag,'_genes_cis_trans_lineages.csv')) %>% as.data.frame()
  dim(cistrans_anno)
  rownames(cistrans_anno) <- NULL
  # cistrans_anno$gene_type <- paste0("in ",cistrans_anno$gene_type)
  summary(as.factor(cistrans_anno$gene_type))
  meta_clone_lg <- data.table::fread(paste0(input_dir,datatag,'_meta_clone_lineages.csv')) %>% as.data.frame()
  meta_clone_lg <- meta_clone_lg[!duplicated(meta_clone_lg$lineage),]
  
  # sum(obs_genes_df$ens_gene_id==rownames(exp_mtx))
  obs_cells_df$lineage <- gsub('Lineage ','L',obs_cells_df$lineage)
  obs_cells_df$lineage <- gsub(': ','-',obs_cells_df$lineage)
  table(obs_genes_df$chromatin_st, obs_genes_df$gene_type)
  
  
  ## Adding new labels
  # save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  # meta_gm_SA609 <- data.table::fread(paste0(save_fig_dir,'meta_gene_module_labels_SA609.csv'))
  # # meta_gm_SA609$gm_manuscript_lb
  # meta_gm_SA609
  # obs_genes_df$gene_type_origin[1:3]
  # obs_genes_df <- obs_genes_df %>%
  #   dplyr::left_join(meta_gm_SA609, by=c('gene_type'='gene_type_module')) %>%
  #   dplyr::rename(gene_type_origin=gene_type, gene_type=gm_manuscript_lb)
  # summary(as.factor(obs_genes_df$gene_type))
  # data.table::fwrite(obs_genes_df, paste0(output_dir,"obs_genes_hm.csv.gz"))
  
  save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  chromatin_df <- data.table::fread(paste0(save_fig_dir, "Gene_Module_Chromatin_Status_v2.csv")) %>% 
    as.data.frame() %>%
    dplyr::filter(datatag=='SA609')
  
  chromatin_df$Module <- paste0("Module",chromatin_df$Module)
  unique(obs_genes_df$gene_type)
  unique(chromatin_df$Module)
  chr_st <- chromatin_df$chromatin_status
  names(chr_st) <- chromatin_df$Module
  obs_genes_df$chromatin_st <- chr_st[obs_genes_df$gene_type]
  summary(as.factor(obs_genes_df$gene_type))
  # obs_genes_df <- obs_genes_df %>%
  #   dplyr::left_join(meta_gm_SA609, by=c('gene_type'='gene_type_module')) %>%
  #   dplyr::rename(gene_type_origin=gene_type, gene_type=gm_manuscript_lb)
  summary(as.factor(obs_genes_df$gene_type))
  obs_genes_df <- get_pretty_gene_type_labels(obs_genes_df)
  summary(as.factor(obs_genes_df$gene_type))
  summary(as.factor(obs_genes_df$chromatin_st))
  output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq_v4/'
  # obs_genes_df$gene_type[1]
  dev.off()
  gap <- viz_genes_exp_lineages_cistrans_anno_hm(as.matrix(exp_mtx), obs_genes_df, obs_cells_df, 
                                                 cistrans_anno, meta_clone_lg,
                                                 output_dir, plttitle)
  # gap
  
  results_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/plots/'
  # saveRDS(gap, paste0(results_dir, datatag,'_hm_plt.rds'))
  gap <- readRDS(paste0(results_dir, datatag,'_hm_plt.rds'))
  gap
  # p <- viz_genes_exp_lineages_hm(as.matrix(exp_mtx), obs_genes_df, obs_cells_df, output_dir, plttitle)
  # p
  
  
  # library(ComplexHeatmap)
  # phm_total <- cowplot::plot_grid(hm_plg, NULL, ncol=2, rel_widths = c(1,0.1),labels = c('c','d'))
  # phm_total
  output_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/"
  # gap <- readRDS(paste0(output_dir,'Activated_Repressed_Transient_genes.rds'))
  # gap
  
  dev.off()
  gap_hmSA609 <- grid.grabExpr(ComplexHeatmap::draw(gap, annotation_legend_side = "bottom",
                                             heatmap_legend_side = "bottom",
                                             padding = unit(c(1, 1, 1, 1), "mm")))
  
  # hm_plg_cistrans <- readRDS(paste0(output_dir,"gene_type_hm_legend.rds"))
  
  res <- plot_legends_heatmap(output_dir)
  # res$hm_plg_cistrans
  
  output_dir <- save_figs_dir
  res <- plot_legends_heatmap(output_dir)
  # ggsave(paste0(save_figs_dir,datatag,"_cistrans_legends.svg"),
  #        plot = cowplot::ggdraw() + cowplot::draw_plot(res$hm_plg_cistrans),
  #        height = 2.5,
  #        width = 1,
  #        # useDingbats=F,
  #        dpi=20)
  # 
  # ggsave(paste0(save_figs_dir,datatag,"_chromatin_legends.svg"),
  #        plot = cowplot::ggdraw() + cowplot::draw_plot(res$hm_plg_chromatin),
  #        height = 2,
  #        width = 1.5,
  #        # useDingbats=F,
  #        dpi=20)
  # ggsave(paste0(save_figs_dir,datatag,"_exp_legends.svg"),
  #        plot = cowplot::ggdraw() + cowplot::draw_plot(res$hm_plg_exp),
  #        height = 1,
  #        width = 2.5,
  #        # useDingbats=F,
  #        dpi=20)
  
  # cowplot::ggdraw() + cowplot::draw_plot(res$hm_plg_chromatin)
  ## in case necessary for yaxis label, otherwise, just add a text
  # yaxis_label = ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_text(c('Pseudotime'), 
  #                                                                                 location = 0.2, rot = 0,
  #                                                                                 just = "center",which='column'))
  # hm_yaxis_title <- grid.grabExpr(ComplexHeatmap::draw(yaxis_label))
  # cowplot::plot_grid(hm_yaxis_title, ncol=1)
  
  # gdh1 <- grid.grabExpr(ComplexHeatmap::draw(gdh, annotation_legend_side = "bottom",
  #                                            heatmap_legend_side = "bottom",
  #                                            padding = unit(c(3, 3, 3, 3), "mm")))
  # phm <- cowplot::plot_grid(gap1, gdh1, rel_widths = c(1,1), ncol=2)
  
  # pgene <- cowplot::plot_grid(exg2, exg2, exg2, exg2, exg1, ncol=5, rel_widths = c(1,1,1,1,1.5))
  # pgene <- cowplot::plot_grid(plotlist = plt_genes_demo, ncol=2) #, labels = c('a','b','c','d')
  # pgene_lg <- cowplot::plot_grid(exg_lg, NULL, ncol=2, rel_widths = c(1.4,0.6)) #,labels = c('f','')
  
  
  
  
}


plot_pathways <- function(){
  ## pathway bootstrap analysis
  # source('/home/htran/Projects/farhia_project/drug_resistant_material/scripts/pipeline/utils/pathway_utils.R')
  script_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/scripts/'
  source(paste0(script_dir, "pipeline/utils/pathway_utils.R"))
  
  ## library(dplyr)
  # save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  # datatag <- 'SA609'
  # pathway_stat <- data.table::fread(paste0(save_fig_dir,'pw_stats_modules_3series.csv')) %>% as.data.frame()
  # View(pathway_stat)
  # pathway_stat <- pathway_stat %>%
  #   dplyr::filter(datatag=='SA609') %>%
  #   dplyr::rename(pathway=reference_set, nb_signf_genes=nb_signif_genes)
  # pathway_stat$gene_type_module <- gsub('Module','M',pathway_stat$gene_type_module)
  # my_font <- "Helvetica"
  # pathway_stat$pathway <- tolower(gsub('HALLMARK_','',pathway_stat$pathway))
  # pathway_stat$datatag <- pathway_stat$gene_type_module
  # pw_609 <- viz_pathways(pathway_stat, datatag, save_fig_dir)
  # pw_609
  # save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  # pw_609 <- readRDS(paste0(save_fig_dir, datatag, "_trajectory_pathways.rds"))
  
  # unique(obs_genes_df$gene_type)
  # unique(obs_genes_df$gene_type_origin)
  library(dplyr)
  library(ggplot2)
  
  
  save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  pathway_stat <- data.table::fread(paste0(save_fig_dir,'pw_stats_modules_3series_gprofiler.csv')) %>% as.data.frame()
  pathway_stat <- pathway_stat %>%
    # dplyr::filter(datatag=='SA609') %>%
    dplyr::rename(pathway=reference_set, nb_signf_genes=nb_signif_genes)
  unique(pathway_stat$gene_type_module)
  
  
  my_font <- "Helvetica"
  pathway_stat$pathway <- tolower(gsub('HALLMARK_','',pathway_stat$pathway))
  pathway_stat1 <- pathway_stat %>%
    dplyr::filter(grepl('SA609',datatag) & p_value <0.05)
  
  
  
  ## Adding new labels
  save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  meta_gm_SA609 <- data.table::fread(paste0(save_fig_dir,'meta_gene_module_labels_SA609.csv'))
  # meta_gm_SA609$gm_manuscript_lb
  meta_gm_SA609
  summary(as.factor(pathway_stat1$gene_type_module))
  pathway_stat1 <- pathway_stat1 %>%
    dplyr::left_join(meta_gm_SA609, by=c('gene_type_module')) %>%
    dplyr::rename(gene_type_origin=gene_type_module, gene_type_module=gm_manuscript_lb)
  summary(as.factor(pathway_stat1$gene_type_module))
  head(pathway_stat1)
  dim(pathway_stat1)
  
  pathway_stat1$gene_type_module <- gsub('Module','M',pathway_stat1$gene_type_module)
  ## TO DO: add new gene module names
  p_pathway <- viz_pathways_barplot(pathway_stat1)
  p_pathway
  results_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/plots/'
  saveRDS(p_pathway, paste0(results_dir,datatag, "p_pathway.rds"))
  p_pathway_SA609 <- readRDS(paste0(results_dir,datatag, "p_pathway.rds"))
  dev.off()
}

plot_smooth_expression <- function(){
  script_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/scripts/trajectory_analysis/'
  datatag <- 'SA609'
  source(paste0(script_dir, "tradeseq_utils.R"))
  input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/SA609/'
  output_dir <- paste0(dirname(input_dir),'/figs_v5/')
  res_smooth_exp <- viz_gene_modules_density_plot(xplt='pseudotime',yplt='gene_exp',
                                colorplt='lineage',datatag='SA609',input_dir, output_dir)
  
  res_smooth_exp$p_smooth_lineage_exp
  ggsave(paste0(results_dir,"Fig6_trajectory_SA609_smooth_exp.svg"),
         plot = res_smooth_exp$p_smooth_lineage_exp,
         height = 7,
         width = 2.5,
         # useDingbats=F,
         dpi=150)
  ## Reloading data
  # p <- readRDS(paste0(output_dir,"density_plot_",datatag, ".rds"))
  # plg <- readRDS(paste0(output_dir,"density_plot_",datatag, "_lg.rds"))
  # res_smooth_exp <- list(p_smooth_lineage_exp=p, plg=plg)
  
}

plot_total <- function(){
  # lg1 <- cowplot::plot_grid(plg, hm_plg, rel_heights = c(3,1), nrow=2)
  # lg1 <- cowplot::plot_grid(hm_plg_exp, hm_plg_cistrans, rel_heights = c(1,1), nrow=2)+
  #   theme(plot.background = element_rect(fill = "white", colour = "white"))
  # pgene_lg <- cowplot::plot_grid(lg1, pw_609, rel_widths = c(0.9,3), nrow=1)
  # pgenes_ls <- cowplot::plot_grid(pgene, pgene_lg, ncol=1, rel_heights = c(2,1))
  # 
  # hm_total <- cowplot::plot_grid(NULL, NULL, gap1, nrow=3, rel_heights = c(0.2, 0.5, 10),
  #                                hjust = -0.1,label_fontface="plain",
  #                                labels=c(' ','Heatmap x: Gene Modules/ y: Pseudotime','')) + 
  #   theme(plot.background = element_rect(fill = "white", colour = "white"))
  # phm_total <- cowplot::plot_grid(hm_total, pgenes_ls, ncol=2, rel_widths = c(1,1),labels = c('c','d'))
  # 
  # p609 <- cowplot::plot_grid(ptraj, phm_total, nrow=2, rel_heights = c(1,1.5))
  # png(paste0(results_dir,"Fig7_trajectory_",datatag,".png"), height = 2*1200, width=2*950,res = 2*72)
  # print(p609)
  # dev.off()
  
  ##
  hm_total_609 <- cowplot::plot_grid(NULL, gap_hmSA609, nrow=2, rel_heights = c(0.03, 5),
                                 hjust = -0.1,label_fontface="plain") + #labels=c(' ','Heatmap x: Gene Modules/ y: Pseudotime','')
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  # p_traj_hm <- cowplot::plot_grid(ptotal_traj_609, hm_total, 
  #                            nrow=3, rel_heights = c(1,1.5), labels=c('a','b'))
  
  p_plg <- cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_plot(res$hm_plg_exp), 
                              cowplot::ggdraw() + cowplot::draw_plot(res$hm_plg_cistrans), 
                              cowplot::ggdraw() + cowplot::draw_plot(res$hm_plg_chromatin),
                              ncol=3, rel_widths = c(1.2, 1, 1),
                                      hjust = -0.1,label_fontface="plain") + 
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  # res_smooth_exp$plg, 
  p_pathway_exp_609 <- cowplot::plot_grid(res_smooth_exp$p_smooth_lineage_exp, NULL, ncol=2, rel_widths = c(1, 1.5),
                     hjust = -0.1, labels=c(' ',' '))
  
  p609 <- cowplot::plot_grid(NULL, ptotal_traj_609, hm_total_609, p_pathway_exp_609, ncol=1, 
                             rel_heights = c(0.02,0.99,1.2,1.5),
                     hjust = -0.1, labels=c(' ',' ',' ', ' ')) +
  theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  results_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/plots/'
  ggsave(paste0(results_dir,"Fig6_trajectory_SA609.png"),
         plot = p609,
         height = 14,
         width = 6,
         # useDingbats=F,
         dpi=150)
  
  ggsave(paste0(results_dir,"Fig6_trajectory_SA609_pathways.png"),
         plot = p_pathway,
         height = 5,
         width = 3.5,
         # useDingbats=F,
         dpi=150)
  ggsave(paste0(results_dir,"Fig6_trajectory_SA609_treatment_lg.png"),
         plot = ts_color,
         height = 2.5,
         width = 3,
         # useDingbats=F,
         dpi=150)
  ggsave(paste0(results_dir,"Fig6_trajectory_SA609_plg_lineage.png"),
         plot = res_smooth_exp$plg,
         height = 1,
         width = 4,
         # useDingbats=F,
         dpi=150)
  ggsave(paste0(results_dir,"Fig6_trajectory_SA609_pathways.png"),
         plot = p_pathway_SA609,
         height = 5,
         width = 3.5,
         # useDingbats=F,
         dpi=150)
  
  final_plt <- cowplot::plot_grid(p609,  p535,NULL, nrow=1, rel_widths = c(1,1,0.01))+ 
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  ggsave(paste0(results_dir,"Fig6_trajectory_total.png"),
         plot = final_plt,
         height = 15,
         width = 12,
         # useDingbats=F,
         dpi=150)
  # ggsave(paste0(results_dir,"Fig7_trajectory_SA609.pdf"),
  #        plot = p609,
  #        height = 13.5,
  #        width = 12,
  #        useDingbats=F,
  #        dpi=150)
  
  
  
}


plot_prevalence <- function(){
  # ## Prevalence plot
  # pts <- readRDS(paste0(results_dir,datatag,'_treatment_trajectory_prevalence.rds'))
  # pc <- readRDS(paste0(results_dir,datatag,'_clone_trajectory_prevalence.rds'))
  # # gb_heatmap = grid.grabExpr(draw(Heatmap(mat)))
  # # The value of padding should be a unit vector with length of four. The four values correspond to 
  # # the space at the bottom, left, top and right sides.
  # prevalence_ts <- grid.grabExpr(ComplexHeatmap::draw(pts, annotation_legend_side = "bottom",
  #                                                     heatmap_legend_side = "bottom",
  #                                                     padding = unit(c(8, 2, 2, 2), "mm")))
  # prevalence_clone <- grid.grabExpr(ComplexHeatmap::draw(pc, annotation_legend_side = "bottom",
  #                                                        heatmap_legend_side = "bottom",
  #                                                        padding = unit(c(8, 2, 2, 2), "mm")))
  # pvalenc <- cowplot::plot_grid(prevalence_ts, prevalence_clone, rel_widths = c(1,1.1), ncol=2)
  # ptraj <- cowplot::plot_grid(ptotal, pvalenc,
  #                             rel_widths = c(3,2.2), ncol=2, labels = c('a','b'))
  # png(paste0(results_dir,"trajectory_",datatag,"_part1.png"), height = 2*400, width=2*900,res = 2*72)
  # print(ptraj)
  # dev.off()
  
}



plot_example_gene <- function(){
  ## ## Specific genes in each module plot
  ## Output of function: viz_given_gene_exp_lineages(obs_genes_ls, meta_genes, ts_sce, figs_dir, datatag)
  
  results_dir2 <- paste0(output_dir,'tradeSeq_3000/SA609_tradeSeq/')
  # demo_genes <- c('TOP2A','COX6C','ID1','CCDC26','HIST1H4C')
  demo_genes <- c('TOP2A','COX6C','ID1','CCDC26')
  # titles <- c('TOP2A: Module 1','COX6C: Module 4','ID1: Module 2','CCDC26: Module 5')
  titles <- c('TOP2A: M1','COX6C: M4','ID1: M2','CCDC26: M5')
  
  plt_genes_demo <- list()
  for(i in seq_len(length(demo_genes))){
    g <- demo_genes[i]
    exg <- readRDS(paste0(results_dir,g,'.rds'))
    exg <- exg + theme(plot.title = element_text(color="black", size=13, hjust = 0.5, family=my_font))
    exg <- exg + labs(title = titles[i])
    if(i%%2==0){
      exg <- exg + theme(legend.position = 'none')
    }
    plt_genes_demo[[g]] <- exg
  }
  
  pgene <- cowplot::plot_grid(plotlist = plt_genes_demo, ncol=2) #, labels = c('a','b','c','d')
  pgene
  
}

get_significant_genes_module5 <- function(){
  library(dplyr)
  ## Get the list of important genes
  output_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/'
  obs_genes_df <- data.table::fread(paste0(output_dir,"SA609_total_genes_modules_act_repr_trans_08_Dec.csv.gz")) %>% as.data.frame()
  dim(obs_genes_df)
  summary(as.factor(obs_genes_df$gene_type_module))
  obs_genes_df <- obs_genes_df %>% 
    dplyr::filter(gene_type_module=='Module5')
  
  input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/cancer_reference_genes/'
  cr <- data.table::fread(paste0(input_dir, 'cisplatin_resistance_corrected.csv')) %>% as.data.frame()
  dim(cr)
  cr$gene_symb  
  
  cf <- data.table::fread(paste0(input_dir, 'CoreFitness_corrected.csv')) %>% as.data.frame()
  cf$gene_symb
  
  bs <- data.table::fread(paste0(input_dir, 'BroadSanger.csv')) %>% as.data.frame()
  dim(bs)
  cosmic <- data.table::fread(paste0(input_dir, 'cosmic.csv')) %>% as.data.frame()
  dim(cosmic)
  ref_genes <- union(cr$gene_symb, cf$gene_symb)
  obs_genes_df$gene_symbol[obs_genes_df$gene_symbol %in% cr$gene_symb]
  obs_genes_df$gene_symbol[obs_genes_df$gene_symbol %in% bs$gene_symb]
  obs_genes_df$gene_symbol[obs_genes_df$gene_symbol %in% cosmic$gene_symb]
  gmt_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/h.all.v7.0.symbols.gmt'
  Hs.H <- fgsea::gmtPathways(gmt_fn) 
  Hs.H
  Hs.H$HALLMARK_TNFA_SIGNALING_VIA_NFKB
  ls_detected_genes <- c()
  for(pw in names(Hs.H)){
    
    detected_genes <- obs_genes_df$gene_symbol[obs_genes_df$gene_symbol %in% Hs.H[[pw]]]
    if(length(detected_genes)>0){
      print(pw)
      print(detected_genes)  
      ls_detected_genes <- c(ls_detected_genes, detected_genes)
    }
    
  }
  ls_detected_genes <- unique(ls_detected_genes)
  signif_genes <- data.frame(gene_symb=ls_detected_genes)
  data.table::fwrite(signif_genes, paste0(output_dir, 'hallmark_signif_genes_module5.csv.gz'))

  
}



get_chromatin_status <- function(){
  # Extracting genes to put into manuscript
  library(dplyr)
  input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/'
  series <- c('SA609','SA535','SA1035')
  pts <- c('Pt4','Pt5','Pt6')
  names(pts) <- series
  chromatin <- data.table::fread(paste0(input_dir,'Gene_Module_Chromatin_Status_v2.csv'))
  # head(chromatin)
  View(chromatin)
  chromatin$chromatin_status <- chromatin$`Promoter Status`
  stat <- tibble::tibble()
  for(datatag in series){
    tag <- datatag
    df <- data.table::fread(paste0(input_dir, datatag, '_total_genes_modules_act_repr_trans_08_Dec.csv.gz'))
    print(dim(df))
    # print(colnames(df))
    df <- df %>%
      dplyr::select(ens_gene_id, gene_symbol, gene_type_module, gene_type, chr, description)
    if(datatag %in% c('SA609','SA535')){
      meta_gm <- data.table::fread(paste0(input_dir,'meta_gene_module_labels_',datatag,'.csv'))
      df <- df %>%
        dplyr::left_join(meta_gm, by=c('gene_type_module')) %>%
        dplyr::select(-gene_type_module) %>%
        dplyr::rename(gene_type_module=gm_manuscript_lb)
    }
    chro <- chromatin %>%
      dplyr::filter(datatag==tag) %>%
      dplyr::mutate(gene_type_module=paste0('Module',Module)) %>%
      dplyr::select(gene_type_module, chromatin_status)
    unique(df$gene_type_module)
    df <- df %>%
      dplyr::left_join(chro, by=c('gene_type_module')) 
    df$patient <- datatag
    print(dim(df))
    stat <- dplyr::bind_rows(stat, df)
  }
  
  dim(stat)
  colnames(stat)
  data.table::fwrite(stat, paste0(input_dir, 'pseudotime_genes_modules_Pt4_Pt5_Pt6.csv.gz'))
  
  
}
