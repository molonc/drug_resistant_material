

# Finding the scripts for 3 DLP plots
# Finding scripts for 3 UMAP plots 
# Adding 10x proportions - clone prevalence
# Adding DLP proportions - clone prevalence



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
  
  series <- c('SA501','SA530','SA604')
  pts_lb <- c('Pt1','Pt2','Pt3') 
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





## Plot 10x umap landscape untreated 3 pdxs
##
plot_10x_UMAP <- function()
{  
  
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  source(paste0(base_dir, 'scripts/drug_manuscript/viz_umap_figs/viz_umaps.R'))
  input_dir <- paste0(base_dir, 'materials/umap_figs/')
  output_dir <- paste0(base_dir, 'materials/umap_figs/figs_rna/')
  # dir.create(output_dir)
  
  ##---------------------------------------------------------------------------------------
  ##-----------------------------------SA501 10x UMAP, barplot prevalence clonealign-------
  datatag <- 'SA501'
  basename <- datatag
  clonealign_stat <- data.table::fread(paste0(input_dir,'clonealign_labels.csv.gz'))
  
  # unique(clonealign_stat$unique_clone)
  clonealign_stat <- clonealign_stat %>%
    dplyr::filter(datatag==basename)%>%
    dplyr::select(cell_id, unique_clone)%>%
    dplyr::rename(clone=unique_clone)
  colnames(clonealign_stat)
  # rm(umap_df)
  umap_df <- data.table::fread(paste0(input_dir,datatag,'_norm_umap.csv.gz')) %>% as.data.frame()
  umap_df$clone <- NULL
  dim(umap_df)
  umap_df$cell_id[1]
  umap_df$scell_id <- paste0(umap_df$id,'_',umap_df$Barcode)
  sum(umap_df$cell_id %in% clonealign_stat$cell_id)
  umap_df <- umap_df %>% left_join(clonealign_stat, by=c('cell_id'))
  summary(as.factor(umap_df$clone))
  
  color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v3.csv.gz')
  cols_use <- get_color_clones(datatag, color_fn) # predefined clone color code for DLP and 10x inferred clones clonealign
  
  tps <- unique(umap_df$timepoint)
  print(tps)
  rna_umap_SA501 <- list()
  obs_treatment <- 'UnRx'
  for(obs_passage in tps){
    plottitle= paste0('10x - ',obs_treatment,': ',obs_passage)
    res <- viz_umap_obs_clones(umap_df, cols_use, datatag, 
                               output_dir, obs_treatment, obs_passage, 
                               plottitle, F)
    rna_umap_SA501[[paste0(obs_treatment,'_',obs_passage)]] <- res
  }
  # rna_umap_SA501$UnRx_X2$p
  # dev.off()
  saveRDS(rna_umap_SA501, paste0(output_dir, datatag, "_umap_plt.rds"))
  # res_SA501_UnRx <- readRDS(paste0(input_dir,'figs/',datatag,'_dlp_prevalence.rds'))
  
  # res_dlp_SA501_lg <- readRDS(paste0(input_dir,'figs/',datatag,'_clone_legend_plt.rds'))
  
  ## Clone prevalence of inferred clones - output of clonealign
  res_prop10x_SA501 <- plot_fill_barplot_wholedataset_rnaseq(umap_df, cols_use, output_dir, 
                                                             datatag, plottitle=NULL, plotlegend=F)
  # res_prop10x_SA501$p
  # res_prop10x_SA501$plg
  saveRDS(res_prop10x_SA501, paste0(output_dir, datatag, "_prevalence_clone_clonealign10x.rds"))
  
  ##---------------------------------------------------------------------------------------
  ##-----------------------------------SA530 10x UMAP, barplot prevalence clonealign-------
  input_dir <- paste0(base_dir, 'materials/umap_figs/')
  output_dir <- paste0(base_dir, 'materials/umap_figs/figs_rna/')
  datatag <- 'SA530'
  tag <- 'Pt2'
  basename <- datatag
  clonealign_stat <- data.table::fread(paste0(input_dir,'clonealign_labels.csv.gz'))
  
  # unique(clonealign_stat$unique_clone)
  clonealign_stat <- clonealign_stat %>%
    dplyr::filter(datatag==basename)%>%
    dplyr::select(cell_id, unique_clone)%>%
    dplyr::rename(clone=unique_clone)
  
  umap_df <- data.table::fread(paste0(input_dir,datatag,'_norm_umap.csv.gz')) %>% as.data.frame()
  dim(umap_df)
  # umap_df$scell_id <- paste0(umap_df$id,'_',umap_df$Barcode)
  # sum(umap_df$cell_id %in% clonealign_stat$cell_id)
  # umap_df <- umap_df %>% left_join(clonealign_stat, by=c('cell_id'))
  
  color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v3.csv.gz')
  cols_use <- get_color_clones(datatag, color_fn) # predefined clone color code for DLP and 10x inferred clones clonealign
  
  tps <- unique(umap_df$timepoint)
  print(tps)
  rna_umap_SA530 <- list()
  obs_treatment <- 'UnRx'
  
  for(obs_passage in tps){
    plottitle= paste0('10x - ',obs_treatment,': ',obs_passage)
    res <- viz_umap_obs_clones(umap_df, cols_use, datatag, 
                               output_dir, obs_treatment, obs_passage, 
                               plottitle, F)
    rna_umap_SA530[[paste0(obs_treatment,'_',obs_passage)]] <- res
  }
  # rna_umap_SA530$UnRx_X3$p
  saveRDS(rna_umap_SA530, paste0(output_dir, datatag, "_umap_plt.rds"))
  
    ## Clone prevalence of inferred clones - output of clonealign
  res_prop10x_SA530 <- plot_fill_barplot_wholedataset_rnaseq(umap_df, cols_use, output_dir, 
                                                             datatag, plottitle=NULL, plotlegend=F)
  # res_prop10x_SA530$p
  # res_prop10x_SA530$plg
  # dev.off()
  saveRDS(res_prop10x_SA530, paste0(output_dir, datatag, "_prevalence_clone_clonealign10x.rds"))
  ##---------------------------------------------------------------------------------------
  ##-----------------------------------SA604 10x UMAP, barplot prevalence clonealign-------
  
  datatag <- 'SA604'
  tag <- 'Pt3'
  basename <- datatag
  clonealign_stat <- data.table::fread(paste0(input_dir,'clonealign_labels.csv.gz'))
  
  # unique(clonealign_stat$unique_clone)
  clonealign_stat <- clonealign_stat %>%
    dplyr::filter(datatag==basename)%>%
    dplyr::select(cell_id, unique_clone)%>%
    dplyr::rename(clone=unique_clone)
  colnames(clonealign_stat)
  
  rm(umap_df)
  umap_df <- data.table::fread(paste0(input_dir,datatag,'_norm_umap.csv.gz')) %>% as.data.frame()
  dim(umap_df)
  umap_df$cell_id[1]
  # umap_df$scell_id <- paste0(umap_df$id,'_',umap_df$Barcode)
  # sum(umap_df$cell_id %in% clonealign_stat$cell_id)
  clonealign_stat$cell_id[1]
  # umap_df <- umap_df %>% left_join(clonealign_stat, by=c('cell_id'))
  umap_df$cell_id <- paste0(umap_df$id,'_', umap_df$Barcode)
  sum(umap_df$cell_id %in% clonealign_stat$cell_id)
  umap_df <- umap_df %>% left_join(clonealign_stat, by=c('cell_id'))
  umap_df <- umap_df %>% 
    dplyr::rename(treatmentSt=treat)
  
  color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v3.csv.gz')
  cols_use <- get_color_clones(datatag, color_fn) # predefined clone color code for DLP and 10x inferred clones clonealign
  
  tps <- unique(umap_df$timepoint)
  print(tps)
  rna_umap_SA604 <- list()
  obs_treatment <- 'UnRx'
  for(obs_passage in tps){
    plottitle= paste0('10x - ',obs_treatment,': ',obs_passage)
    res <- viz_umap_obs_clones(umap_df, cols_use, datatag, 
                               output_dir, obs_treatment, obs_passage, 
                               plottitle, F)
    rna_umap_SA604[[paste0(obs_treatment,'_',obs_passage)]] <- res
  }
  # rna_umap_SA604$UnRx_X6$p # please exclude X6 from analysis
  saveRDS(rna_umap_SA604, paste0(output_dir, datatag, "_umap_plt.rds"))
  ## Clone prevalence of inferred clones - output of clonealign
  res_prop10x_SA604 <- plot_fill_barplot_wholedataset_rnaseq(umap_df, cols_use, output_dir, 
                                                             datatag, plottitle=NULL, plotlegend=F)
  
  # dev.off()
  saveRDS(res_prop10x_SA604, paste0(output_dir, datatag, "_prevalence_clone_clonealign10x.rds"))
  
  
  
  
}


plot_DLP_barplot <- function(){
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  source(paste0(base_dir, 'scripts/drug_manuscript/viz_umap_figs/viz_umaps.R'))
  input_dir <- paste0(base_dir, 'materials/cell_clones/')
  output_dir <- paste0(base_dir, 'materials/umap_figs/figs_rna/')
  
  ##---------------------------------------------------------------------------------------
  ##-----------------------------------SA501 barplot prevalence from sitka tree -----------
  
  datatag <- 'SA501'
  library_grouping_fn <- paste0(input_dir, datatag, '_DLP_library_groupings.csv.gz')
  cell_clones_fn <- paste0(input_dir, datatag, '_cell_clones.csv.gz')
  cell_clones <- data.table::fread(cell_clones_fn) %>% as.data.frame()
  dim(cell_clones)
  # cols_use <- make_clone_palette(unique(cell_clones$clone_id)) # old version
  color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v3.csv.gz')
  cols_use <- get_color_clones(datatag, color_fn) # predefined clone color code for DLP and 10x inferred clones clonealign
  obs_treatment <- 'UnRx'
  obs_passage <- 'X2'
  cell_clones$timepoint <- 'X2'
  
  cell_clones <- get_meta_data(cell_clones_fn, library_grouping_fn, datatag)
  cell_clones$treatment_desc <- paste0('DLP-',cell_clones$treatment_desc)
  res_barDLP_501 <- plot_fill_barplot_wholedataset_v2(cell_clones, cols_use, output_dir, 
                                                      datatag, plottitle=NULL, plotlegend=F)
  
  # res_barDLP_501$plg
  # res_barDLP_501$p
  # dev.off()
  saveRDS(res_barDLP_501, paste0(output_dir, datatag, "_barplot_DLP.rds"))
  
  
  ##---------------------------------------------------------------------------------------
  ##-----------------------------------SA530 barplot prevalence from sitka tree -----------
  datatag <- 'SA530'
  library_grouping_fn <- paste0(input_dir, datatag, '_DLP_library_groupings.csv.gz')
  cell_clones_fn <- paste0(input_dir, datatag, '_cell_clones.csv.gz')
  cell_clones <- data.table::fread(cell_clones_fn) %>% as.data.frame()
  dim(cell_clones)
  # cols_use <- make_clone_palette(unique(cell_clones$clone_id)) # old version
  color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v3.csv.gz')
  cols_use <- get_color_clones(datatag, color_fn) # predefined clone color code for DLP and 10x inferred clones clonealign
  obs_treatment <- 'UnRx'
  obs_passage <- 'X3'
  cell_clones$timepoint <- 'X3'
  cell_clones <- get_meta_data(cell_clones_fn, library_grouping_fn, datatag)
  cell_clones$treatment_desc <- paste0('DLP-',cell_clones$treatment_desc)
  res_barDLP_530 <- plot_fill_barplot_wholedataset_v2(cell_clones, cols_use, output_dir, 
                                                      datatag, plottitle=NULL, plotlegend=F)
  
  # res_barDLP_530$p
  # res_barDLP_530$plg
  saveRDS(res_barDLP_530, paste0(output_dir, datatag, "_barplot_DLP.rds"))
  
  
  
  
  ##---------------------------------------------------------------------------------------
  ##-----------------------------------SA604 barplot prevalence from sitka tree -----------
  input_dir <- paste0(base_dir, 'materials/cell_clones/')
  output_dir <- paste0(base_dir, 'materials/umap_figs/figs_rna/')
  
  datatag <- 'SA604'
  library_grouping_fn <- paste0(input_dir, datatag, '_DLP_library_groupings.csv.gz')
  cell_clones_fn <- paste0(input_dir, datatag, '_cell_clones.csv.gz')
  cell_clones <- data.table::fread(cell_clones_fn) %>% as.data.frame()
  dim(cell_clones)
  # cols_use <- make_clone_palette(unique(cell_clones$clone_id)) # old version
  color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v3.csv.gz')
  cols_use <- get_color_clones(datatag, color_fn) # predefined clone color code for DLP and 10x inferred clones clonealign
  cell_clones <- get_meta_data(cell_clones_fn, library_grouping_fn, datatag)
  cell_clones$treatment_desc <- paste0('DLP-',cell_clones$treatment_desc)
  res_barDLP_604 <- plot_fill_barplot_wholedataset_v2(cell_clones, cols_use, output_dir, 
                                                      datatag, plottitle=NULL, plotlegend=F)
  # res_barDLP_604$p
  # res_barDLP_604$plg

  saveRDS(res_barDLP_604, paste0(output_dir, datatag, "_barplot_DLP.rds"))
  
}


plot_SUPP_fig1 <- function(){
  
  ## Eric
  ## Change base_dir to your github folder, ex: yourdir/ddrug_resistant_material/materials/
  ## Loading rds files and do plotting
  ## You can use plot_grid or save each patient data into 1 image and paste into svg file, up to you
  ## 
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  source(paste0(base_dir, 'scripts/drug_manuscript/cnv_viz_utils.R'))
  input_dir <- paste0(base_dir, 'materials/umap_figs/figs_rna/')

  
  ## Loading DLP plots from folder or run functions above
  series <- c('SA501','SA530','SA604')
  pts_lb <- c('Pt1','Pt2','Pt3') 
  names(pts_lb) <- series
  res_501 <- readRDS(paste0(input_dir, series[1],'_DLP_heatmap.rds')) # contain main plots and legends
  res_530 <- readRDS(paste0(input_dir, series[2],'_DLP_heatmap.rds'))
  res_604 <- readRDS(paste0(input_dir, series[3],'_DLP_heatmap.rds'))
  
  ## if you have dlp_heatmap_ls from above script, or you can reload it here. 
  dlp_heatmap_ls <- list()
  dlp_heatmap_ls[['SA501']] <- res_501
  dlp_heatmap_ls[['SA530']] <- res_530
  dlp_heatmap_ls[['SA604']] <- res_604
  
  
  supp_fig1_leftside <- cowplot::plot_grid(dlp_heatmap_ls$SA501$cnv_plot, dlp_heatmap_ls$SA530$cnv_plot, 
                                  dlp_heatmap_ls$SA604$cnv_plot, res_501$plg,
                                  ncol = 1, rel_heights = c(1,1,1, 0.2)) + #align = 'v', labels = c('SA501','SA530','SA604','SA609','SA535','SA1035')
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  
  # supp_fig1_leftside
  
  ggsave(paste0(save_dir,"SUPPFig1_3untreated_series_medianCNV.pdf"),
         plot = supp_fig1_leftside,
         height = 10,
         width = 6,
         useDingbats=F, # from Sam's suggestion
         dpi = 150
  )
  
  ggsave(paste0(save_dir,"SUPPFig1_3untreated_series.png"),
         plot = supp_fig1_leftside,
         height = 10,
         width = 6,
         type = "cairo-png",
         dpi=150
  )
  
  
  
  ## SUPP Fig 1, Eric
  res_barDLP_501 <- readRDS(paste0(output_dir, series[1], "_barplot_DLP.rds"))
  res_barDLP_530 <- readRDS(paste0(output_dir, series[2], "_barplot_DLP.rds"))
  res_barDLP_604 <- readRDS(paste0(output_dir, series[3], "_barplot_DLP.rds"))
  
  res_prop10x_SA501 <- readRDS(paste0(output_dir, series[1], "_prevalence_clone_clonealign10x.rds"))
  res_prop10x_SA530 <- readRDS(paste0(output_dir, series[2], "_prevalence_clone_clonealign10x.rds"))
  res_prop10x_SA604 <- readRDS(paste0(output_dir, series[3], "_prevalence_clone_clonealign10x.rds"))
  
  
  rna_umap_SA501 <- readRDS(paste0(output_dir, series[1], "_umap_plt.rds"))
  rna_umap_SA530<- readRDS(paste0(output_dir, series[2], "_umap_plt.rds"))
  rna_umap_SA604 <- readRDS(paste0(output_dir, series[3], "_umap_plt.rds"))
  
  clone_plg_501 <- plot_clone_color_legend('SA501', base_dir, ncols_grid=1)
  p501_total <- cowplot::plot_grid(res_barDLP_501$p, res_prop10x_SA501$p, clone_plg_501, rna_umap_SA501$UnRx_X2$p, 
                                   nrow=1, rel_widths = c(1,1,1,2)) + 
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  clone_plg_530 <- plot_clone_color_legend('SA530', base_dir, ncols_grid=1)
  p530_total <- cowplot::plot_grid(res_barDLP_530$p, res_prop10x_SA530$p, clone_plg_530, rna_umap_SA530$UnRx_X3$p, 
                                   nrow=1, rel_widths = c(1,1,1,2)) + 
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  # p530_total
  
  clone_plg_604 <- plot_clone_color_legend('SA604', base_dir, ncols_grid=3)
  ##rna_umap_SA604$UnRx_X6 ## noted exclude X6 - 10x data from analysis
  
  p604_part1 <- cowplot::plot_grid(res_barDLP_604$p, res_prop10x_SA604$p, clone_plg_604, rel_heights = c(1,1, 1), ncol=1)
  p604_part2 <- cowplot::plot_grid(rna_umap_SA604$UnRx_X7$p, rna_umap_SA604$UnRx_X8$p, 
                                   rna_umap_SA604$UnRx_X9$p, ncol=2) + 
    theme(plot.background = element_rect(fill = "white", colour = "white")) 
  p604_total <- cowplot::plot_grid(p604_part1, p604_part2, rel_widths = c(1,1))+ 
    theme(plot.background = element_rect(fill = "white", colour = "white"))

    
  supp_fig1_rightside <- cowplot::plot_grid(p501_total,p530_total, p604_total, rel_heights = c(1,1,3), ncol=1)+ #labels=pts_lb, 
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  ## The best way is defining a layout using plot_grid function but it takes time...
  
  supp_fig1 <- cowplot::plot_grid(supp_fig1_leftside, supp_fig1_rightside, rel_widths = c(1,1), ncol=2)+ 
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  ggsave(paste0(output_dir,"SUPPFig1_3untreated_series.pdf"),
         plot = supp_fig1,
         height = 11,
         width = 13,
         useDingbats=F, # from Sam's suggestion
         dpi = 150
  )
  
  ggsave(paste0(output_dir,"SUPPFig1_3untreated_series.png"),
         plot = supp_fig1,
         height = 11,
         width = 13,
         type = "cairo-png",
         dpi=150
  )
  
  
} 


