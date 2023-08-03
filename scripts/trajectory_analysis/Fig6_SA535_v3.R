script_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/scripts/'
# source(paste0(script_dir, "trajectory_analysis/slingshot_utils.R"))
source(paste0(script_dir, "trajectory_analysis/tradeseq_utils.R"))
library(grid)

plot_heatmap_gene_modules <- function(){
  
  ## Heatmap average expression plot
  ## Extracting heatmap of average genes expression x lineages 
  datatag <- 'SA535'
  save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  genes_df <- data.table::fread(paste0(save_dir, 'SA535_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
  dim(genes_df)
  colnames(genes_df)
  # summary(as.factor(genes_df$gene_type_module))
  genes_df <- genes_df %>%
    dplyr::select(-gene_type) %>%
    dplyr::rename(gene_type=gene_type_module)
  summary(as.factor(genes_df$gene_type))
  plttitle <- 'Activated Repressed Transient genes'
  save_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/"
  # ts_sce <- readRDS(paste0(save_dir,'tradeSeq_3000/', "fitGAM_out.rds"))
  # dim(ts_sce)
  # phm <- viz_heatmap(ts_sce, genes_df, save_dir, datatag, plttitle)
  output_dir <- save_dir
  
  output_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/SA535/'
  exp_mtx <- data.table::fread(paste0(output_dir,"mtx_hm.csv.gz")) %>% as.data.frame()
  obs_genes_df <- data.table::fread(paste0(output_dir,"obs_genes_hm.csv.gz")) %>% as.data.frame()
  obs_cells_df <- data.table::fread(paste0(output_dir,"obs_cells_hm.csv.gz")) %>% as.data.frame()
  rownames(exp_mtx) <- exp_mtx$ens_gene_id
  exp_mtx$ens_gene_id <- NULL
  dim(exp_mtx)
  dim(obs_genes_df)
  dim(obs_cells_df)
  input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/cis_trans_lineages/')
  cistrans_anno <- data.table::fread(paste0(input_dir, datatag,'_genes_cis_trans_lineages.csv')) %>% as.data.frame()
  dim(cistrans_anno)
  # head(cistrans_anno)
  rownames(cistrans_anno) <- NULL
  # cistrans_anno$gene_type <- paste0("in ",cistrans_anno$gene_type)
  meta_clone_lg <- data.table::fread(paste0(input_dir,datatag,'_meta_clone_lineages.csv')) %>% as.data.frame()
  meta_clone_lg <- meta_clone_lg[!duplicated(meta_clone_lg$lineage),]
  meta_clone_lg$lineage_desc <- gsub('Lineage ','L',meta_clone_lg$lineage_desc)
  meta_clone_lg$lineage_desc <- gsub(': ','-',meta_clone_lg$lineage_desc)
  meta_clone_lg$lineage_desc <- gsub(': ','-',meta_clone_lg$lineage_desc)
  unique(meta_clone_lg$lineage_desc)
  meta_clone_lg$lineage_desc <- ifelse(meta_clone_lg$lineage_desc=="L4-Rx,1Rx","L4-1Rx",meta_clone_lg$lineage_desc)
  # sum(obs_genes_df$ens_gene_id==rownames(exp_mtx))
  unique(obs_genes_df$gene_type)
  # head(obs_cells_df)
  obs_cells_df$lineage <- gsub(': ','-',obs_cells_df$lineage)
  unique(obs_cells_df$lineage)
  
  save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  chromatin_df <- data.table::fread(paste0(save_fig_dir, "Gene_Module_Chromatin_Status_v2.csv")) %>% 
    as.data.frame() %>%
    dplyr::filter(datatag=='SA535')
  
  chromatin_df$Module <- paste0("Module",chromatin_df$Module)
  unique(obs_genes_df$gene_type)
  unique(chromatin_df$Module)
  chr_st <- chromatin_df$chromatin_status
  names(chr_st) <- chromatin_df$Module
  obs_genes_df$chromatin_st <- chr_st[obs_genes_df$gene_type]
  summary(as.factor(obs_genes_df$chromatin_st))
  
  summary(as.factor(obs_genes_df$gene_type))
  obs_genes_df <- get_pretty_gene_type_labels(obs_genes_df)
  summary(as.factor(obs_genes_df$gene_type))
  obs_cells_df$lineage <- ifelse(obs_cells_df$lineage=="L4-Rx,1Rx","L4-1Rx",obs_cells_df$lineage)
  gap <- viz_genes_exp_lineages_cistrans_anno_hm(as.matrix(exp_mtx), obs_genes_df, obs_cells_df,
                                                 cistrans_anno, meta_clone_lg,
                                                 output_dir, plttitle)
  output_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/"
  # gap <- readRDS(paste0(output_dir,'Activated_Repressed_Transient_genes.rds'))
  
  results_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/plots/'
  # saveRDS(gap, paste0(results_dir, datatag,'_hm_plt.rds'))
  gap <- readRDS(paste0(results_dir, datatag,'_hm_plt.rds'))
  gap
  dev.off()
  gap_hmSA535 <- grid.grabExpr(ComplexHeatmap::draw(gap, annotation_legend_side = "bottom",
                                             heatmap_legend_side = "bottom",
                                             padding = unit(c(1, 1, 1, 1), "mm")))
  
  # hm_plg_cistrans <- readRDS(paste0(output_dir,"gene_type_hm_legend.rds"))
  
  res <- plot_legends_heatmap(output_dir)
  hm_total_535 <- cowplot::plot_grid(NULL, gap_hmSA535, nrow=2, rel_heights = c(0.1, 5),
                                 hjust = -0.1,label_fontface="plain") + #labels=c(' ','Heatmap x: Gene Modules/ y: Pseudotime','')
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  
}

plot_trajectory <- function(){
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
                   panel.spacing = unit(c(0, 0, 0, 0), "null")) + 
    labs(title='L4-1Rx')
  
  l4 <- readRDS(paste0(results_dir,'ts_slingshot_out_Lineage4_',datatag,'.rds'))
  l4 <- l4 + theme(plot.title = element_text(color="black", size=11, hjust = 0.5, family=my_font),
                   plot.margin = unit(c(0, 0, 0, 0), "null"),
                   panel.spacing = unit(c(0, 0, 0, 0), "null"))
  
  total <- readRDS(paste0(results_dir,'ts_slingshot_out_wholedataset_',datatag,'.rds'))
  total <- total + theme(plot.title = element_text(color="black", size=16, hjust = 0, family=my_font, face="bold"))
                         #axis.title = element_text(color="black", size=10, hjust = 0.5, family=my_font))
  # "/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/figs_v4/"
  ts_color <- readRDS(paste0(results_dir,'treatment_plt.rds'))
  # ts_color <- plot_legend_treatmentSt(results_dir)
  ts_color <- readRDS(paste0(results_dir,'Treatment_Cycles_plt.rds'))
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
  # dev.off()
  
  # TO DO
  ## Adding legends here
  
  # ggsave(paste0(save_dir,"pseudo_gene_modules_",datatag, ".png"),
  #        plot = ptotal_traj_SA535,
  #        height = 4,
  #        width = 5,
  #        # useDingbats=F,
  #        dpi=150)
  # 
  # ggsave(paste0(save_dir,"pseudo_gene_modules_",datatag, "_color_lg.png"),
  #        plot = ts_color,
  #        height = 2.5,
  #        width = 3.5,
  #        # useDingbats=F,
  #        dpi=150)
  
}



plot_smooth_expression <- function(){
    script_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/scripts/trajectory_analysis/'
    source(paste0(script_dir, "tradeseq_utils.R"))
    input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/SA535/'
    output_dir <- paste0(dirname(input_dir),'/figs_v5/')
    
    res_smooth_exp_535 <- viz_gene_modules_density_plot(xplt='pseudotime',yplt='gene_exp',
                                                    colorplt='lineage',datatag='SA535',input_dir, output_dir)
    # p <- readRDS(paste0(output_dir,"density_plot_",datatag, ".rds"))
    # plg <- readRDS(paste0(output_dir,"density_plot_",datatag, "_lg.rds"))
    # res_smooth_exp <- list(p_smooth_lineage_exp=p, plg=plg)
    res_smooth_exp_535$p_smooth_lineage_exp
    res_smooth_exp_535$plg
    ggsave(paste0(output_dir,"Fig6_trajectory_SA535_smooth_exp.svg"),
           plot = res_smooth_exp_535$p_smooth_lineage_exp,
           height = 7,
           width = 2.5,
           # useDingbats=F,
           dpi=150)
    # res_smooth_exp_535 <- list(p_smooth_lineage_exp=p, plg=plg)
}  
plot_pathways <- function(){
  ## pathway bootstrap analysis
  source('/home/htran/Projects/farhia_project/drug_resistant_material/scripts/pipeline/utils/pathway_utils.R')
  library(dplyr)
  library(ggplot2)
  tag <- 'SA535'
  save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  pathway_stat <- data.table::fread(paste0(save_fig_dir,'pw_stats_modules_3series_gprofiler.csv')) %>% as.data.frame()
  pathway_stat <- pathway_stat %>%
    # dplyr::filter(datatag=='SA609') %>%
    dplyr::rename(pathway=reference_set, nb_signf_genes=nb_signif_genes)
  
  my_font <- "Helvetica"
  pathway_stat$pathway <- tolower(gsub('HALLMARK_','',pathway_stat$pathway))
  pathway_stat1 <- pathway_stat %>%
    dplyr::filter(grepl(tag,datatag) & p_value <0.05)
  
  head(pathway_stat1)
  dim(pathway_stat1)
  
  ## Get pretty labels
  save_fig_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/'
  meta_gm <- data.table::fread(paste0(save_fig_dir,'meta_gene_module_labels_',datatag,'.csv'))
  meta_gm
  summary(as.factor(pathway_stat1$gene_type_module))
  pathway_stat1 <- pathway_stat1 %>%
    dplyr::left_join(meta_gm, by=c('gene_type_module')) %>%
    dplyr::rename(gene_type_origin=gene_type_module, gene_type_module=gm_manuscript_lb)
  summary(as.factor(pathway_stat1$gene_type_module))
  # head(pathway_stat1)
  dim(pathway_stat1)
  pathway_stat1$gene_type_module <- gsub('Module','M',pathway_stat1$gene_type_module)
  ## TO DO: add new gene module names
  p_pathway_SA535 <- viz_pathways_barplot(pathway_stat1)
  p_pathway_SA535
  dev.off()
}

plot_total <- function(){
  p_plg <- cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_plot(res$hm_plg_exp), 
                              cowplot::ggdraw() + cowplot::draw_plot(res$hm_plg_cistrans), 
                              cowplot::ggdraw() + cowplot::draw_plot(res$hm_plg_chromatin),
                              ncol=3, rel_widths = c(1.2, 1, 1),
                              hjust = -0.1,label_fontface="plain")
  # res_smooth_exp$plg, 
  p_pathway_exp_535 <- cowplot::plot_grid(res_smooth_exp_535$p_smooth_lineage_exp, NULL, ncol=2, rel_widths = c(1, 1.5),
                                      hjust = -0.1, labels=c(' ',''))
  
  p535 <- cowplot::plot_grid(NULL, ptotal_traj_SA535, hm_total_535, p_pathway_exp_535, nrow=4, 
                             rel_heights = c(0.01,1,1.15,1.45),
                             hjust = -0.1, labels=c(' ',' ',' ', '')) + #, labels=c(NULL,'b','d', NULL)
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  results_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/plots/'
  ggsave(paste0(results_dir,"Fig6_trajectory_SA535_pathways.png"),
         plot = p_pathway_SA535,
         height = 5,
         width = 3.5,
         # useDingbats=F,
         dpi=150)
  ggsave(paste0(results_dir,"Fig6_trajectory_SA535_plg_lineage.png"),
         plot = res_smooth_exp_535$plg,
         height = 1,
         width = 6,
         # useDingbats=F,
         dpi=150)
  
  ggsave(paste0(results_dir,"Fig6_trajectory_SA535.png"),
         plot = p535,
         height = 14,
         width = 7,
         # useDingbats=F,
         dpi=150)
  
  
  ggsave(paste0(results_dir,"Fig6_trajectory_SA535.pdf"),
         plot = p535,
         height = 13.5,
         width = 12,
         useDingbats=F,
         dpi=150)
  ggsave(paste0(results_dir,"Fig6_trajectory_SA535_treatment_lg.png"),
         plot = ts_color,
         height = 2.5,
         width = 3,
         # useDingbats=F,
         dpi=150)
}  



save_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/figs/'
ggsave(paste0(save_dir,"umaps_lineages_SA535.png"),
       plot = p,
       height = 4,
       width = 6.5,
       # useDingbats=F,
       dpi=150)
save_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/figs/'
ggsave(paste0(save_dir,"umaps_lineages_plg_SA535.png"),
       plot = lg2,
       height = 2,
       width = 4.5,
       # useDingbats=F,
       dpi=150)
