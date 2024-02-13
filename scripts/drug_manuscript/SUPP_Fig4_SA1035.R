

##---------------------------------------------------------------------------------------
##-------------------------------median CNV profiles - from sitka output-------
plot_medianCNV_heatmap <- function(){
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  source(paste0(base_dir, 'scripts/drug_manuscript/cnv_viz_utils.R'))
  input_dir <- paste0(base_dir, 'materials/cell_clones/')
  save_dir <- paste0(input_dir, 'fig_medianCNV/')
  output_dir <- paste0(base_dir, 'materials/umap_figs/figs_rna/')
  if(!dir.exists(save_dir)){
    dir.create(save_dir) 
  }
  
  series <- c('SA1035')
  pts_lb <- c('Pt6') 
  names(pts_lb) <- series
  
  dlp_heatmap_ls <- list()
  for(datatag in series){
    copynumber_fn <- paste0(input_dir, datatag, '_total_merged_filtered_states.csv.gz')
    cellclone_fn <- paste0(input_dir, datatag, '_cell_clones.csv.gz')
    library_grouping_fn <- paste0(input_dir, datatag, '_DLP_library_groupings.csv')
    df_cnv <- get_median_genotype_v3(copynumber_fn, datatag, save_dir,
                                     cellclone_fn, library_grouping_fn) 
    print(dim(df_cnv))
    res <- plot_CNV_profile(df_cnv, clones= levels(df_cnv$clone),
                            plttitle=paste0(pts_lb[datatag],' DLP - median copy number profile'))
    dlp_heatmap_ls[[datatag]] <- res
    saveRDS(res, paste0(output_dir, datatag, '_DLP_heatmap.rds'))
    # res$cnv_plot # main plot
    # res$plg # legend
    # dev.off()
  }
  
  
  
}

plot_coefficients_fitness <- function(){
  
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  source(paste0(base_dir, 'scripts/drug_manuscript/coefficient_viz_utils.R'))
  output_dir <- paste0(base_dir, 'materials/umap_figs/figs_rna/')
  # color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v3.csv.gz')
  color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v4.csv.gz')
  coef_df <- data.table::fread(paste0(base_dir,'materials/fitness_paper_DLP/s_coeff_filter.csv.gz'))
  
  unique(coef_df$datatag)
  ######### Untreated ######### 
  dat <- coef_df %>%
    dplyr::filter(datatag=='SA1035U')                  
  dim(dat)
  
  tag <- 'SA1035'
  plottitle <- 'UnRx'
  p_coeff_UnRx <- plot_box_posterior_v2(dat, color_fn, tag, plottitle, ymax=1.4, ymin=0.6,
                                        clone_axis_dir = 'bottom', sort_by = 'median', add_linebreak = FALSE)
  
  # p_coeff_UnRx
  saveRDS(p_coeff_UnRx, paste0(output_dir, tag, "_UnRx_coeff_fitness.rds"))
  
  ######### Treated ######### 
  dat <- coef_df %>%
    dplyr::filter(datatag=='SA1035T')                  
  dim(dat)
  
  tag <- 'SA1035'
  plottitle <- 'Rx'
  p_coeff_Rx <- plot_box_posterior_v2(dat, color_fn, tag, plottitle, ymax=1.4, ymin=0.6, 
                                      clone_axis_dir = 'bottom', sort_by = 'median', add_linebreak = FALSE)
  
  saveRDS(p_coeff_Rx, paste0(output_dir, tag, "_Rx_coeff_fitness.rds"))
  
}
plot_10x_UMAP <- function()
{  
  
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  source(paste0(base_dir, 'scripts/drug_manuscript/viz_umap_figs/viz_umaps.R'))
  
  input_dir <- paste0(base_dir, 'materials/umap_figs/')
  output_dir <- paste0(base_dir, 'materials/umap_figs/figs_rna/')
  # output_dir <- paste0(base_dir, 'materials/umap_figs/testing/')
  # dir.create(output_dir)
  datatag <- 'SA1035'
  basename <- datatag
  
  ## Revision, change clonealign results
  clonealign_dir <- paste0('/home/htran/Projects/farhia_project/drug_resistant_material/materials/clonealign_plot/clonealign/SA1035-v6/')
  clonealign_stat <- data.table::fread(paste0(clonealign_dir,'SA1035-v6_clonealign_labels_total.csv.gz'))
  dim(clonealign_stat)
  ## First submission
  # t <- data.table::fread(paste0(input_dir,'clonealign_labels.csv.gz'))
  # dim(t)
  # unique(t$clone)
  # t <- t %>%
  #   dplyr::filter(datatag=="SA1035" & clone!="unassigned")
  # dim(clonealign_stat)
  # unique(clonealign_stat$clone)
  # clonealign_stat <- clonealign_stat %>%
  #   dplyr::filter(clone!="unassigned")
  
  # sum(t$cell_id %in% umap_df$cell_id)
  clonealign_stat <- data.table::fread(paste0(input_dir,'clonealign_labels.csv.gz'))
  # unique(clonealign_stat$unique_clone)
  clonealign_stat <- clonealign_stat %>%
    dplyr::filter(datatag==basename)%>%
    dplyr::select(cell_id, unique_clone)%>%
    dplyr::rename(clone=unique_clone)
  # colnames(clonealign_stat)
  # dim(clonealign_stat)
  # clonealign_stat$cell_id[1]
  # summary(as.factor(clonealign_stat$clone))
  # rm(umap_df)

  # umap_df <- data.table::fread(paste0(input_dir,datatag,'_norm_umap.csv.gz')) %>% as.data.frame()
  umap_df <- data.table::fread(paste0(input_dir,datatag,'_norm_umap_filtered_outliers.csv.gz')) %>% as.data.frame()
  dim(umap_df)
  summary(as.factor(umap_df$clone))
  umap_df$clone <- NULL ## revision, update new data
  umap_df$cell_id[1]
  clonealign_stat$cell_id[1]
  # sum(clonealign_stat$cell_id %in% umap_df$cell_id)
  dim(clonealign_stat)
  umap_df <- umap_df %>%
    dplyr::left_join(clonealign_stat, by='cell_id') %>%
    dplyr::mutate(clone = 
                    case_when(
                      is.na(clone) ~ 'None', 
                      TRUE ~ clone
                    ))
  summary(as.factor(umap_df$treatmentSt))
  summary(as.factor(umap_df$clone))
  sids <- sapply(strsplit(umap_df$cell_id,'_'), function(x){
    return(x[1])
  })  
  umap_df$sample_id <- as.character(sids)
  # table(umap_df$treatmentSt, umap_df$sample_id)
  # table(umap_df$clone, umap_df$sample_id)
  # unique(umap_df$sample_id)
  # umap_df <- umap_df %>%
  #   dplyr::filter(timepoint=='X4')
  
  # umap_df <- umap_df %>%
  #   dplyr::filter(timepoint=='X4')
  # dim(umap_df)
  # summary(as.factor(umap_df$clone))
  umap_df <- umap_df %>%
    dplyr::mutate(clone=get_unique_clone_id(clone))
  
  # summary(as.factor(umap_df$treatmentSt))
  # sum(umap_df$cell_id %in% clonealign_stat$cell_id)
  ## For pseudotime clonal labels, revision manuscript
  # data.table::fwrite(umap_df, '/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/clone_labels_unique_SA1035.csv.gz')
  
  # filtered_cells_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')
  # sce <- readRDS(paste0(output_dir, datatag,'_3000_rd_sce_v2.rds'))
  
  # umaps <- reducedDims(sce)[['UMAP']] %>% as.data.frame()
  # umaps$cluster_label <- sce$cluster_label
  # umaps$cluster_label <- as.factor(umaps$cluster_label)
  # umaps$clone <- sce$clone
  # umaps$treatmentSt <- colData(sce)[rownames(umaps),'treatmentSt']
  # umaps$cell_id <- rownames(umaps)
  # umaps$timepoint <- sce$timepoint
  # unique(umaps$treatmentSt)
  # data.table::fwrite(umaps, paste0(input_dir,datatag,'_norm_umap_filtered_outliers.csv.gz'))
  
  # umap_df$scell_id <- paste0(umap_df$id,'_',umap_df$Barcode)
  # sum(umap_df$cell_id %in% clonealign_stat$cell_id)
  # umap_df <- umap_df %>% left_join(clonealign_stat, by=c('cell_id'))
  
  # color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v3.csv.gz')
  color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v4.csv.gz')
  cols_use <- get_color_clones(datatag, color_fn) # predefined clone color code for DLP and 10x inferred clones clonealign
  
  
  # umap_df$scell_id <- paste0(umap_df$id,'_',umap_df$Barcode)
  # sum(umap_df$scell_id %in% clonealign_stat$cell_id)
  # umap_df <- umap_df %>% left_join(clonealign_stat, by=c('scell_id'='cell_id'))
  # summary(as.factor(umap_df$clone))
  res_prop10x_SA1035 <- plot_fill_barplot_wholedataset_rnaseq(umap_df, cols_use, output_dir, 
                                                              datatag, plottitle=NULL, plotlegend=F, facet_order='vertical')
  
  # dev.off()
  # res_prop10x_SA1035$p
  saveRDS(res_prop10x_SA1035, paste0(output_dir, datatag, "_prevalence_clone_clonealign10x.rds"))
  
  # dim(umap_df)
  tps <- gtools::mixedsort(unique(umap_df$timepoint))
  unique(tps)
  rna_SA1035 <- list()
  # obs_treatment <- 'UnRx'
  for(obs_treatment in c('UnRx','Rx','RxH')){
    # res <- viz_umap_obs_clones(umap_df, cols_use, datatag, 
    #                            output_dir, obs_treatment, NULL, 
    #                            obs_treatment, F)
    # rna_SA535[[paste0(obs_treatment,'_total')]] <- res
    for(obs_passage in tps){
      plottitle= paste0(obs_treatment,': ',obs_passage)
      res <- viz_umap_obs_clones(umap_df, cols_use, datatag, 
                                 output_dir, obs_treatment, obs_passage, 
                                 plottitle, F)
      
      rna_SA1035[[paste0(obs_treatment,'_',obs_passage)]] <- res
    }
  }
  # res$p
  saveRDS(rna_SA1035, paste0(output_dir,datatag,"_umap_plt.rds"))
  # rna_SA1035 <- readRDS(paste0(output_dir,datatag,"_umap_plt.rds"))
  # rna_SA1035$UnRx_X5$p
  # clone_plg <- plot_clone_color_legend(datatag, base_dir, ncols_grid=2)
  
  # p1035_total <- cowplot::plot_grid(rna_SA1035$UnRx_X5$p,rna_SA1035$UnRx_X6$p, rna_SA1035$UnRx_X7$p, rna_SA1035$UnRx_X8$p,
                                    # rna_SA1035$Rx_X5$p, rna_SA1035$Rx_X6$p,rna_SA1035$Rx_X7$p, rna_SA1035$Rx_X8$p,
                                    # clone_plg, rna_SA1035$RxH_X6$p, rna_SA1035$RxH_X7$p, rna_SA1035$RxH_X8$p,
                                    # ncol=4, align = 'hv')
  
  # unique(umap_df$treatmentSt)
  
  
}


plot_DLP_barplot <- function(){
  ##---------------------------------------------------------------------------------------
  ##-----------------------------------SA604 barplot prevalence from sitka tree -----------
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  input_dir <- paste0(base_dir, 'materials/cell_clones/')
  output_dir <- paste0(base_dir, 'materials/umap_figs/figs_rna/')
  
  datatag <- 'SA1035'
  library_grouping_fn <- paste0(input_dir, datatag, '_DLP_library_groupings.csv.gz')
  cell_clones_fn <- paste0(input_dir, datatag, '_cell_clones.csv.gz')
  cell_clones <- get_meta_data(cell_clones_fn, library_grouping_fn, datatag)
  unique(cell_clones$treatment_desc)
  
  summary(as.factor(cell_clones$clone_id))
  # cell_clones <- data.table::fread(cell_clones_fn) %>% as.data.frame()
  # dim(cell_clones)
  # obs_samples <- c('SA609X3XB01584','SA609X4XB03080','SA609X5XB03223',
  #                  'SA609X6XB03447','SA609X7XB03554','SA609X4XB003083',
  #                  'SA609X5XB03230','SA609X6XB03404','SA609X7XB03505',
  #                  'SA609X5XB03231','SA609X6XB03401','SA609X7XB03510')
  # cell_clones <- cell_clones %>%
  #   dplyr::filter(sample_id %in% obs_samples)
  # cell_clones$sample_id
  # # 
  # cell_clones_fn <- paste0(input_dir, datatag, '_cell_clones_main_line.csv.gz')
  # data.table::fwrite(cell_clones, cell_clones_fn)
  # unique(cell_clones$sample)
  # dim(cell_clones)
  # unique(cell_clones$treatment_desc)
  # cols_use <- make_clone_palette(unique(cell_clones$clone_id)) # old version
  # color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v3.csv.gz')
  color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v4.csv.gz')
  cols_use <- get_color_clones(datatag, color_fn) # predefined clone color code for DLP and 10x inferred clones clonealign
  
  
  cell_clones$treatment_desc <- paste0('DLP - ',cell_clones$treatment_desc)
  res_barDLP_1035 <- plot_fill_barplot_wholedataset_v2(cell_clones, cols_use, output_dir, 
                                                      datatag, plottitle=NULL, plotlegend=F, facet_order='vertical')
  # res_barDLP_1035$p
  # res_barDLP_1035$plg
  
  saveRDS(res_barDLP_1035, paste0(output_dir, datatag, "_barplot_DLP.rds"))
  
  
}




plot_SUPP_fig2 <- function(){
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  input_dir <- paste0(base_dir, 'materials/cell_clones/')
  output_dir <- paste0(base_dir, 'materials/umap_figs/figs_rna/')
  source(paste0(base_dir, 'scripts/drug_manuscript/viz_umap_figs/viz_umaps.R'))
  
  datatag <- 'SA1035'
  
  medianCNV_plt_topplot <- readRDS(paste0(output_dir, datatag, '_DLP_heatmap.rds'))
  # medianCNV_plt_topplot1 <- cowplot::plot_grid(medianCNV_plt_topplot$cnv_plot, medianCNV_plt_topplot$plg, rel_heights = c(1,0.1), ncol=1)
  # top_plt <- cowplot::plot_grid(medianCNV_plt_topplot1, NULL, rel_widths = c(4,2), nrow=1)+ 
  #   theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  
  medianCNV_plt_topplot1 <- cowplot::plot_grid(medianCNV_plt_topplot$cnv_plot, medianCNV_plt_topplot$plg, 
                                               rel_heights = c(1,0.1), ncol=1,
                                               labels=c('a ',NULL),label_size=11, hjust=0)
  
  p_coeff_UnRx <- readRDS(paste0(output_dir, datatag, "_UnRx_coeff_fitness.rds"))
  p_coeff_Rx <- readRDS(paste0(output_dir, datatag, "_Rx_coeff_fitness.rds"))
  plt_topplot2 <-cowplot::plot_grid(NULL, p_coeff_UnRx, p_coeff_Rx,ncol=1, 
                                    rel_heights = c(0.11, 1, 1), labels = c('b Fitness coefficients',NULL,NULL),
                                    label_size=11, hjust=0)
  
  top_plt <- cowplot::plot_grid(medianCNV_plt_topplot1, plt_topplot2, NULL, rel_widths = c(3.8,2.2, 0.1), nrow=1)#+ 
    # theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  
  
  clone_plg <- plot_clone_color_legend(datatag, base_dir, ncols_grid=2)
  
  res_barDLP_1035 <- readRDS(paste0(output_dir, datatag, "_barplot_DLP.rds"))
  
  res_prop10x_SA1035 <- readRDS(paste0(output_dir, datatag, "_prevalence_clone_clonealign10x.rds"))
  
  rna_SA1035 <- readRDS(paste0(output_dir,datatag,"_umap_plt.rds"))
  
  
  p_leftside <- cowplot::plot_grid(rna_SA1035$UnRx_X5$p,rna_SA1035$UnRx_X6$p, rna_SA1035$UnRx_X7$p, rna_SA1035$UnRx_X8$p,
                                   rna_SA1035$Rx_X5$p, rna_SA1035$Rx_X6$p,rna_SA1035$Rx_X7$p, rna_SA1035$Rx_X8$p,
                                   clone_plg, rna_SA1035$RxH_X6$p, rna_SA1035$RxH_X7$p, rna_SA1035$RxH_X8$p,
                                   ncol=4, align = 'hv')#+
    # theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  # p_prevalence <- cowplot::plot_grid(res_prop10x_SA1035$p, res_barDLP_1035$p, rel_widths = c(1, 1), nrow=1)
  p_prevalence <- cowplot::plot_grid(res_prop10x_SA1035$p, NULL, rel_widths = c(1, 1), nrow=1)
  
  p_prevalence1 <- cowplot::plot_grid(NULL, p_prevalence, rel_heights = c(0.04, 1), ncol=1)
  
  p_bottom <- cowplot::plot_grid(p_leftside, NULL, p_prevalence1, rel_widths = c(4,0.05,2), nrow=1) #+ 
  
    # theme(plot.background = element_rect(fill = "white", colour = "white"))
  # p_total <- cowplot::plot_grid(top_plt, p_bottom, rel_heights = c(2,3.5), ncol=1) + 
  #   theme(plot.background = element_rect(fill = "white", colour = "white"))
  # 
  p_bottom1 <- cowplot::plot_grid(NULL, NULL, NULL, labels=c('  c 10x UMAP across treatments & time points',
                                                             '   d 10x inferred \nclones prevalence','    e DLP clones\n prevalence'),
                                  label_size=11, hjust=0, rel_widths = c(4,1,1),nrow=1)#, fontfamily = "Helvetica"
  
  p_total <- cowplot::plot_grid(top_plt, p_bottom1, p_bottom, rel_heights = c(2,0.13,3.5), ncol=1) + 
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  p_total <- cowplot::plot_grid(NULL, NULL, p_bottom, rel_heights = c(2,0.13,3.5), ncol=1) #+ 
    # theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  ggsave(paste0(output_dir,"SUPPFig4_",datatag, "_revision.svg"),
         plot = p_total,
         height = 12,
         width = 11,
         # type = "cairo-png",
         dpi=200
  )
  
  ggsave(paste0(output_dir,"SUPPFig4_",datatag, ".png"),
         plot = p_total,
         height = 12,
         width = 11,
         type = "cairo-png",
         dpi=150
  )
  
}
