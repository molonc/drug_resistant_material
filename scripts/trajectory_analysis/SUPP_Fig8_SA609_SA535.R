
script_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/scripts/trajectory_analysis/'
source(paste0(script_dir, "slingshot_utils.R"))
source(paste0(script_dir, "tradeseq_utils.R"))

## Plot proportion of cells across lineage
## Characterize each lineage based on cells prevalence
characterizing_lineages_SA609 <- function(){
  
  datatag <- 'SA609'
  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
  input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/')
  
  # Load data for trajectory plotting
  output_dir <- save_dir
  nfeatures_use <- 3000
  sce <- readRDS(paste0(save_dir, datatag,'_',nfeatures_use,'_rd_sce.rds'))
  dim(sce)
  sce$sample[1]
  metacell <- colData(sce) %>% as.data.frame()
  data.table::fwrite(metacell, paste0(output_dir, datatag,'_meta_cells.csv.gz'))
  dim(metacell)
  # metacell <- data.table::fread(paste0(output_dir, datatag,'_meta_cells.csv.gz'))
  # sce <- readRDS(paste0(output_dir, datatag,'_',nfeatures_use,'_rd_sce.rds'))
  # pca_df <- data.table::fread(paste0(output_dir, datatag,'_',nfeatures_use, "_norm_pca.csv")) %>% as.data.frame()
  # umap_df <- data.table::fread(paste0(output_dir, datatag,'_',nfeatures_use, "_norm_umap.csv")) %>% as.data.frame()
  # dim(umap_df)
  start_cls <- '10'
  rd_use <- 'PCA'
  nfeatures_use <- 3000
  datatag <- 'SA609'
  crv_umap_embed <- readRDS(paste0(save_dir, "slingshot_",datatag,'_',paste(start_cls, collapse='_'),"_UMAP_embed_crv.rds"))
  input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/cis_trans_lineages/')
  meta_clone_lg <- data.table::fread(paste0(input_dir,datatag,'_meta_clone_lineages.csv')) %>% as.data.frame()
  meta_clone_lg <- meta_clone_lg[!duplicated(meta_clone_lg$lineage),]
  meta_lineage <- meta_clone_lg %>%
    dplyr::select(lineage, lineage_desc)
  p <- plot_trajectory_clones_prevalence_SA609(metacell, input_dir, save_dir, meta_lineage)
  tag <- 'clone'
  p_clone_prevalence_609 <- readRDS(paste0(paste0(output_dir,'figs_v3/'),datatag, "_",tag,"_trajectory_prevalence.rds"))
  tag <- 'treatment'
  p_treatment_prevalence_609 <- readRDS(paste0(paste0(output_dir,'figs_v3/'),datatag, "_",tag,"_trajectory_prevalence.rds"))
  
  
  
}

plot_example_genes_SA609 <- function(){
  library(tidyverse)
  library(ComplexHeatmap)
  obs_genes_ls <- c('TOP2A','COX6C','ID1','CCDC26')
  datatag <- 'SA609'
  figs_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/figs_v3/"
  save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/'
  total_genes <- data.table::fread(paste0(save_dir, 'SA609_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
  dim(total_genes)
  meta_genes <- total_genes %>%
    dplyr::filter(gene_symbol %in% obs_genes_ls) %>%
    dplyr::select(-description)
  
  
  meta_genes <- meta_genes %>% 
    remove_rownames %>%
    dplyr::select(ens_gene_id, gene_symbol) %>%
    tibble::column_to_rownames('gene_symbol')
  
  head(meta_genes)
  datatag <- 'SA609'
  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
  input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/')
  nfeatures_use <- 3000
  output_dir <- paste0(save_dir,'tradeSeq_3000/')
  ts_sce <- readRDS(paste0(save_dir,'tradeSeq_3000/', "fitGAM_out.rds"))
  output_fn <- paste0(save_dir,'tradeSeq_3000/','sigf_gene_exp.csv.gz')
  
  avg_gene_exp <- get_average_gene_exp_per_lineage(datatag, ts_sce, meta_genes$ens_gene_id, output_fn, nEstimatedPoints=100, save_data=T)
  # save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/')
  # ts_sce <- readRDS(paste0(save_dir,'tradeSeq_3000/', "fitGAM_out.rds"))
  dim(avg_gene_exp)
  colnames(avg_gene_exp)
  viz_given_gene_exp_lineages(obs_genes_ls, meta_genes, avg_gene_exp, figs_dir, datatag)
  
  input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'  
  meta_gm <- data.table::fread(paste0(input_dir,'meta_gene_module_labels_',datatag,'.csv.gz'))
  meta_genes <- meta_genes %>% left_join(meta_gm, by='gene_type_module')
  meta_genes$gm_manuscript_lb <- gsub('Module','M',meta_genes$gm_manuscript_lb)
  meta_genes$gm_manuscript_lb <- paste0(meta_genes$gm_manuscript_lb,': ',meta_genes$gene_symbol)
  obs_genes_titles <- meta_genes$gm_manuscript_lb
  names(obs_genes_titles) <- meta_genes$gene_symbol
  obs_genes_titles
  obs_genes_ls <- c('TOP2A','ID1','COX6C','CCDC26')
  
  
  plt_609_ls <- list()
  for(gsymb in obs_genes_ls){
    pg <- readRDS(paste0(figs_dir,gsymb,'.rds'))  
    pg <- pg + labs(title = obs_genes_titles[gsymb])
    plt_609_ls[[gsymb]] <- pg + theme(legend.position = 'none',
                                  plot.title = element_text(size=12, family = my_font, color='black'))
  }  
  pg <- readRDS(paste0(figs_dir,gsymb,'.rds'))  
  pg <- pg + theme(legend.position = 'bottom')
  plg_lg <- cowplot::get_legend(pg)
  plg_609 <- cowplot::ggdraw(plg_lg)
  
  p_obs_genes_609 <- cowplot::plot_grid(plotlist = plt_609_ls, nrow = 2)
  p_obs_genes_609_total <- cowplot::plot_grid(p_obs_genes_609, plg_609, nrow=2, rel_heights=c(2, 0.2))
  # gb <- grid.grabExpr(draw(p1, annotation_legend_side = "bottom",
  #                          heatmap_legend_side = "bottom"),
  #                     padding = unit(c(1, 1, 2, 2), "mm"))#,merge_legend = TRUE
  # 
  # ptotal <- cowplot::plot_grid(gb, p_obs_genes, nrow = 2, rel_heights = c(2,1))
  # png(paste0(output_dir,"reverse_33genes_summary_",datatag,".png"), height = 2*650, width=2*400,res = 2*72)
  # print(ptotal)
  # dev.off()
  # 
  # png(paste0(output_dir,"reverse_33genes_summary_",datatag,".png"), height = 2*550, width=2*650,res = 2*72)
  # print(ptotal)
  # dev.off()
  
}

characterizing_lineages_SA535 <- function(){
  datatag <- 'SA535'
  input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')
  output_dir <- save_dir
  
  start_cls <- '7'
  # start_cls <- '12'
  rd_use <- 'PCA'
  nfeatures_use <- 3000
  crv_umap_embed <- readRDS(paste0(save_dir, "slingshot_",datatag,start_cls,"_UMAP_embed_crv.rds"))
  # crv1 <- readRDS(paste0(save_dir, "slingshot_pseudotime_",datatag,"_",start_cls,'_',rd_use,"_crv.rds"))
  # sce <- readRDS(paste0(output_dir, datatag,'_',nfeatures_use,'_rd_sce.rds'))
  sce <- readRDS(paste0(output_dir, datatag,'_',nfeatures_use,'_rd_sce_clones.rds'))
  dim(sce)
  colnames(sce)[1]
  colnames(sce) <- sce$cell_id
  sce$treatmentSt <- sce$treat
  rownames(sce)[1]
  metacell <- colData(sce) %>% as.data.frame()
  data.table::fwrite(metacell, paste0(output_dir, datatag,'_meta_cells.csv.gz'))
  dim(metacell)
  # metacell <- data.table::fread(paste0(output_dir, datatag,'_meta_cells.csv.gz'))
  
  input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/cis_trans_lineages/')
  meta_clone_lg <- data.table::fread(paste0(input_dir,datatag,'_meta_clone_lineages.csv')) %>% as.data.frame()
  meta_clone_lg <- meta_clone_lg[!duplicated(meta_clone_lg$lineage),]
  meta_clone_lg$lineage_desc <- gsub('Lineage ','L',meta_clone_lg$lineage_desc)
  meta_clone_lg$lineage_desc <- gsub(': ','-',meta_clone_lg$lineage_desc)
  meta_clone_lg$lineage_desc <- gsub(': ','-',meta_clone_lg$lineage_desc)
  unique(meta_clone_lg$lineage_desc)
  meta_clone_lg$lineage_desc <- ifelse(meta_clone_lg$lineage_desc=="L4-Rx,1Rx","L4-1Rx",meta_clone_lg$lineage_desc)
  meta_lineage <- meta_clone_lg %>%
    dplyr::select(lineage, lineage_desc)
  # meta_lineage
  p <- plot_trajectory_clones_prevalence_SA535(metacell, input_dir, output_dir, meta_lineage)
  tag <- 'treatment_clone'
  p_clone_prevalence_535 <- readRDS(paste0(output_dir,'figs_v3/',datatag, "_",tag,"_trajectory_prevalence.rds"))
  tag <- 'treatment'
  p_treatment_prevalence_535 <- readRDS(paste0(output_dir,'figs_v3/',datatag, "_",tag,"_trajectory_prevalence.rds"))
  # p_treatment_prevalence_535
}  
  
plot_example_genes_SA535 <- function(){
  library(tidyverse)
  library(ComplexHeatmap)
  my_font <- "Helvetica"
  ## Visualize specific genes expr
  datatag <- 'SA535'
  figs_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/figs_v3/"
  # 
  save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  total_genes <- data.table::fread(paste0(save_dir, datatag,'_total_genes_modules_act_repr_trans_08_Dec.csv.gz')) %>% as.data.frame()
  dim(total_genes)
  obs_genes_ls <- c('SAA1','CENPW','SMC4','RRM2')
  # obs_genes_ls <- total_genes$ens_gene_id
  meta_genes <- total_genes %>%
    dplyr::filter(gene_symbol %in% obs_genes_ls) %>%
    dplyr::select(-description)
  meta_genes <- meta_genes %>% 
    remove_rownames %>%
    dplyr::select(ens_gene_id, gene_symbol) %>%
    tibble::column_to_rownames('gene_symbol')
  
  input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/'
  ts_sce <- readRDS(paste0(input_dir,'tradeSeq_3000/', "fitGAM_out.rds"))
  dim(ts_sce)
  output_fn <- paste0(input_dir,'tradeSeq_3000/','sigf_gene_exp.csv.gz')
  avg_gene_exp <- get_average_gene_exp_per_lineage(datatag, ts_sce, meta_genes$ens_gene_id, output_fn, nEstimatedPoints=100, save_data=T)
  # save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/')
  obs_genes_ls <- 
  meta_genes <- total_genes %>%
    dplyr::filter(gene_symbol %in% obs_genes_ls)
  dim(avg_gene_exp)
  
  meta_genes <- meta_genes %>% 
    remove_rownames %>%
    dplyr::select(ens_gene_id, gene_symbol) %>%
    tibble::column_to_rownames('gene_symbol')
  unique(avg_gene_exp$lineage_desc)
  # dim(avg_gene_exp)
  # head(meta_genes)
  viz_given_gene_exp_lineages(obs_genes_ls, meta_genes, avg_gene_exp, figs_dir, datatag)
  
  input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'  
  meta_gm <- data.table::fread(paste0(input_dir,'meta_gene_module_labels_',datatag,'.csv.gz'))
  meta_genes <- meta_genes %>% left_join(meta_gm, by='gene_type_module')
  meta_genes$gm_manuscript_lb <- gsub('Module','M',meta_genes$gm_manuscript_lb)
  meta_genes$gm_manuscript_lb <- paste0(meta_genes$gm_manuscript_lb,': ',meta_genes$gene_symbol)
  obs_genes_titles <- meta_genes$gm_manuscript_lb
  names(obs_genes_titles) <- meta_genes$gene_symbol
  obs_genes_ls <- c('SMC4','RRM2','SAA1','CENPW')
  plt_535_ls <- list()
  for(gsymb in obs_genes_ls){
    pg <- readRDS(paste0(figs_dir,gsymb,'.rds'))  
    pg <- pg + labs(title = obs_genes_titles[gsymb])
    pg <- pg + theme(legend.position = 'none',
                     plot.title = element_text(size=12, family = my_font, color='black'))
    plt_535_ls[[gsymb]] <- pg
  }  
  pg <- readRDS(paste0(figs_dir,gsymb,'.rds'))  
  pg <- pg + theme(legend.position = 'bottom')
  plg_lg <- cowplot::get_legend(pg)
  plg_535 <- cowplot::ggdraw(plg_lg)
  p_obs_genes_535 <- cowplot::plot_grid(plotlist = plt_535_ls, nrow = 2)
  p_obs_genes_535_total <- cowplot::plot_grid(p_obs_genes_535, plg_535, nrow=2, rel_heights=c(2, 0.2))
  
}

viz_SUPPFig8 <- function(){
  library(grid)
  datatag <- 'SA609'
  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
  input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/')
  output_dir <- save_dir
  figs_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/figs_v3/"
  ## Prevalence plt
  tag <- 'clone'
  p_clone_prevalence_609 <- readRDS(paste0(paste0(output_dir,'figs_v3/'),datatag, "_",tag,"_trajectory_prevalence.rds"))
  tag <- 'treatment'
  p_treatment_prevalence_609 <- readRDS(paste0(paste0(output_dir,'figs_v3/'),datatag, "_",tag,"_trajectory_prevalence.rds"))
  
  plt_clone_prevalence_609 <- grid::grid.grabExpr(ComplexHeatmap::draw(p_clone_prevalence_609, annotation_legend_side = "bottom",
                           heatmap_legend_side = "right"),
                      padding = unit(c(10, 1, 1, 20), "mm"))#,merge_legend = TRUE
  plt_treatment_prevalence_609 <- grid::grid.grabExpr(ComplexHeatmap::draw(p_treatment_prevalence_609, annotation_legend_side = "bottom",
                                           heatmap_legend_side = "right"),
                      padding = unit(c(10, 1, 1, 20), "mm"))#,merge_legend = TRUE
  
  prevalence_total_609 <- cowplot::plot_grid(NULL, plt_treatment_prevalence_609, plt_clone_prevalence_609, 
                                             rel_widths = c(0.05,1,1.1), nrow = 1)
  # prevalence_total_609
  # cowplot::ggdraw(plt_clone_prevalence_609)
  # p_clone_prevalence_609
  # p_treatment_prevalence_609
  ## Example genes
  # obs_genes_ls <- c('TOP2A','COX6C','ID1','CCDC26')
  # plt_609_ls <- list()
  # for(gsymb in obs_genes_ls){
  #   pg <- readRDS(paste0(figs_dir,gsymb,'.rds'))  
  #   plt_609_ls[[gsymb]] <- pg + theme(legend.position = 'none')
  # }  
  # pg <- readRDS(paste0(figs_dir,gsymb,'.rds')) 
  # pg <- pg + theme(legend.position = 'bottom')
  # plg_lg <- cowplot::get_legend(pg)
  # plg_609 <- cowplot::ggdraw(plg_lg)
  # p_obs_genes_609 <- cowplot::plot_grid(plotlist = plt_609_ls, nrow = 2)
  # p_obs_genes_609_total <- cowplot::plot_grid(p_obs_genes_609, plg_609, nrow=2, rel_heights=c(2, 0.2))
  # 
  # # p_obs_genes_609_total
  # dev.off()
  
  datatag <- 'SA535'
  input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')
  output_dir <- save_dir
  figs_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/figs_v3/"
  
  ## Prevalence plt
  tag <- 'treatment_clone'
  p_clone_prevalence_535 <- readRDS(paste0(output_dir,'figs_v3/',datatag, "_",tag,"_trajectory_prevalence.rds"))
  tag <- 'treatment'
  p_treatment_prevalence_535 <- readRDS(paste0(output_dir,'figs_v3/',datatag, "_",tag,"_trajectory_prevalence.rds"))
  # p_treatment_prevalence_535
  plt_clone_prevalence_535 <- grid.grabExpr(ComplexHeatmap::draw(p_clone_prevalence_535, annotation_legend_side = "bottom",
                                                                 heatmap_legend_side = "right"),
                                            padding = unit(c(10, 1, 1, 20), "mm"))
  plt_treatment_prevalence_535 <- grid.grabExpr(ComplexHeatmap::draw(p_treatment_prevalence_535, annotation_legend_side = "bottom",
                                                                     heatmap_legend_side = "right"),
                                                padding = unit(c(10, 1, 1, 20), "mm"))#,merge_legend = TRUE
  # cowplot::plot_grid(plt_clone_prevalence_535)
  prevalence_total_535 <- cowplot::plot_grid(NULL, plt_treatment_prevalence_535, plt_clone_prevalence_535, 
                                             rel_widths = c(0.05, 1, 1.1), nrow = 1)
  # prevalence_total_535
  
  # p_clone_prevalence_535
  # p_treatment_prevalence_535
  
  ## Example genes
  # obs_genes_ls <- c('SAA1','CENPW','SMC4','RRM2')
  # plt_535_ls <- list()
  # for(gsymb in obs_genes_ls){
  #   pg <- readRDS(paste0(figs_dir,gsymb,'.rds'))  
  #   plt_535_ls[[gsymb]] <- pg + theme(legend.position = 'none')
  # }  
  # pg <- readRDS(paste0(figs_dir,gsymb,'.rds'))  
  # pg <- pg + theme(legend.position = 'bottom')
  # plg_lg <- cowplot::get_legend(pg)
  # plg_535 <- cowplot::ggdraw(plg_lg)
  # p_obs_genes_535 <- cowplot::plot_grid(plotlist = plt_535_ls, nrow = 2)
  # p_obs_genes_535_total <- cowplot::plot_grid(p_obs_genes_535, plg_535, nrow=2, rel_heights=c(2, 0.2))
  # 
  # # p_obs_genes_535_total
  # dev.off()
  
  
  # prevalence_summary <- cowplot::plot_grid(prevalence_total_609, NULL,
  #                                          prevalence_total_535, NULL, rel_widths = c(1,0.05,1,0.05), nrow=1)
  # obs_genes_summary <- cowplot::plot_grid(p_obs_genes_609_total, NULL, p_obs_genes_535_total, NULL, 
  #                                         rel_widths = c(1,0.05,1,0.05), nrow=1)
  
  p609 <- cowplot::plot_grid(prevalence_total_609, NULL,
                             p_obs_genes_609_total, NULL, rel_widths = c(1,0.05,0.8,0.01), nrow=1)
  p535 <- cowplot::plot_grid(prevalence_total_535, NULL,
                             p_obs_genes_535_total, NULL, rel_widths = c(1,0.05,0.8,0.01), nrow=1)
  
  
  # suppFig8 <- cowplot::plot_grid(NULL, prevalence_summary, obs_genes_summary, rel_heights = c(0.05, 1.3,1), nrow=3)
  p_top_legend <- cowplot::plot_grid(NULL, NULL, labels=c(' a Pt4 lineages - treatment      Lineages - treatment clone',
                                                          ' b Pt4 Example gene from gene modules'),
                                  label_size=11, hjust=0, rel_widths = c(1.05,0.81),nrow=1)#, fontfamily = "Helvetica"
  p_middle_legend <- cowplot::plot_grid(NULL, NULL, labels=c(' c Pt5 lineages - treatment      Lineages - treatment clone',
                                                          ' d Pt5 Example gene from gene modules'),
                                     label_size=11, hjust=0, rel_widths = c(1.05,0.81),nrow=1)#, fontfamily = "Helvetica"
  
  suppFig8 <- cowplot::plot_grid(p_top_legend, p609, p_middle_legend, p535, NULL, rel_heights = c(0.07, 1, 0.15, 1.2, 0.15), nrow=5)
  
  
  suppFig8_part1 <- cowplot::plot_grid(prevalence_total_609, NULL, prevalence_total_535, ncol=1, rel_heights = c(1, 0.15, 1.2))
  
  
  save_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/trajectory_genes/figs/'
  ggsave(paste0(save_dir,"SUPP_Fig8_SA609_SA535_trajectory.svg"),
         plot = suppFig8,
         height = 13,
         width = 11,
         # useDingbats=F,
         dpi=250)
  
  ggsave(paste0(save_dir,"SUPP_Fig8_SA609_SA535_trajectory_part1.svg"),
         plot = suppFig8_part1,
         height = 12,
         width = 7,
         # useDingbats=F,
         dpi=250)
  
  
  
}


checking_pathway_Farhia <- function(){
  library(dplyr)
  save_fig_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'
  tag <- 'SA609'
  pathway_stat <- data.table::fread(paste0(save_fig_dir,'pw_stats_modules_3series_gprofiler.csv')) %>% as.data.frame()
  pathway_stat <- pathway_stat %>%
    dplyr::filter(datatag==tag)
  View(pathway_stat)
}
## Add some example of genes with different lineage