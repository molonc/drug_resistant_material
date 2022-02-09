plot_SA609_same_chips <- function(output_dir){
  input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/'
  output_dir <- paste0(input_dir,'figs/')
  datatag <- 'SA609'
  pca_df <- data.table::fread(paste0(input_dir,datatag,'_norm_pca.csv')) %>% as.data.frame()
  head(pca_df)
  pca_df$PC_1
  pca_df$clone <- factor(pca_df$clone, levels = sort(unique(pca_df$clone)))
  p <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color=clone)) + 
    geom_point(size=1.1, alpha=0.8) +
    # geom_point(shape = 1,size = 1.5, colour = "grey")
    scale_color_manual(values=cols_use, name='') + 
    labs(title = plottitle) +
    xlim(xl[1], xl[2]) + 
    ylim(yl[1], yl[2])
  
  p <- p + thesis_theme
  
  
  datatag <- 'SA609_samechip'
  umap_df <- data.table::fread(paste0(input_dir,datatag,'_norm_umap.csv')) %>% as.data.frame()
  head(umap_df)
  cols_use <- make_clone_palette(unique(umap_df$clone))
  unique(umap_df$timepoint)
  obs_treatment <- 'Rx'
  obs_passage <- 'X5'
  plottitle= paste0(obs_treatment,': ',obs_passage)
  res1 <- viz_umap_obs_clones(umap_df, cols_use, datatag, output_dir, obs_treatment, obs_passage, plottitle, F)
  
  obs_treatment <- 'Rx'
  obs_passage <- 'X6'
  plottitle= paste0(obs_treatment,': ',obs_passage)
  res2 <- viz_umap_obs_clones(umap_df, cols_use, datatag, output_dir, obs_treatment, obs_passage, plottitle, F)
  
  obs_treatment <- 'Rx'
  obs_passage <- 'X7'
  plottitle= paste0(obs_treatment,': ',obs_passage)
  res3 <- viz_umap_obs_clones(umap_df, cols_use, datatag, output_dir, obs_treatment, obs_passage, plottitle, F)
  
  p609_samechip <- cowplot::plot_grid(res1$p,res2$p, res3$p, nrow=1)
  p609_total <- cowplot::plot_grid(p609_samechip, res1$plg, ncol=1, rel_heights = c(1, 0.1))
  
  png(paste0(output_dir,datatag,".png"), height = 2*250, width=2*800,res = 2*72)
  print(p609_total)
  dev.off()
  
}  
get_clone_labels <- function(output_dir){
  # Get clone labels first
  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/'
  input_dir <- '/home/htran/Projects/farhia_project/drug_resistance/differential_expression/comps/input_data/'
  clonealign_dir <- paste0(base_dir,'clonealign/')
  dlp_dir <- paste0(base_dir,'dlp/')
  fns <- list.files(clonealign_dir)    
  fns <- fns[grepl('*csv',fns)]
  
  clonealign_stat <- tibble::tibble()
  for(f in fns){
    # First clonealign output
    sdf <- data.table::fread(paste0(clonealign_dir,f)) %>% as.data.frame()
    # print(dim(sdf))
    clonealign_stat <- dplyr::bind_rows(clonealign_stat, sdf)
  }  
  dim(clonealign_stat)
  clonealign_stat$library_id <- gsub('.cache/','',clonealign_stat$Sample[1])
  clonealign_stat$library_id <- gsub('/filtered_feature_bc_matrix','',clonealign_stat$library_id)
  clonealign_stat$Sample <- NULL
  clonealign_stat$datatag <- stringr::str_sub(clonealign_stat$id,1,6)
  clonealign_stat$datatag <- gsub('X$','',clonealign_stat$datatag)
  unique(clonealign_stat$datatag)
  clonealign_stat$lcell_id <- paste0(clonealign_stat$library_id,'_',clonealign_stat$Barcode)
  clonealign_stat$unique_clone <- get_unique_clone_id(clonealign_stat$clone)
  output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/'
  data.table::fwrite(clonealign_stat, paste0(output_dir,'clonealign_labels.csv'))
  # View(head(clonealign_stat))
}  


plot_SA535 <- function(output_dir){
  output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/figs_rna/'
  # dir.create(output_dir)
  input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/'
 
  datatag <- 'SA535'
  clonealign_stat <- data.table::fread(paste0(dirname(output_dir),'/clonealign_labels.csv'))
  
  unique(clonealign_stat$unique_clone)
  clonealign_stat <- clonealign_stat %>%
    dplyr::filter(datatag=='SA535')%>%
    dplyr::select(cell_id, unique_clone)%>%
    dplyr::rename(clone=unique_clone)
  
  
  umap_df <- data.table::fread(paste0(input_dir,datatag,'_norm_umap.csv')) %>% as.data.frame()
  head(umap_df[1:3,1:3])
  umap_df$cell_id[1]
  clonealign_stat$lcell_id[1]
  dim(umap_df)
  umap_df$scell_id <- paste0(umap_df$id,'_',umap_df$Barcode)
  sum(umap_df$scell_id %in% clonealign_stat$cell_id)
  umap_df <- umap_df %>% left_join(clonealign_stat, by=c('scell_id'='cell_id'))
  summary(as.factor(umap_df$clone))
  
  cols_use <- get_color_clones(datatag)
  tps <- unique(umap_df$timepoint)
  rna_SA535 <- list()
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
      rna_SA535[[paste0(obs_treatment,'_',obs_passage)]] <- res
    }
  }
  
  #rna_SA535$UnRx_X5$p,
  p535_total <- cowplot::plot_grid(rna_SA535$UnRx_X6$p, rna_SA535$UnRx_X7$p,rna_SA535$UnRx_X8$p, rna_SA535$UnRx_X9$p,NULL,res_SA535_UnRx$p,
                                   rna_SA535$Rx_X6$p,rna_SA535$Rx_X7$p, rna_SA535$Rx_X8$p,rna_SA535$Rx_X9$p,rna_SA535$Rx_X10$p,res_SA535_Rx$p,
                                   res_dlp_SA535_lg$plg, rna_SA535$RxH_X7$p, rna_SA535$RxH_X8$p,rna_SA535$RxH_X9$p,rna_SA535$RxH_X10$p,res_SA535_RxH$p,
                                   ncol=6, align = 'hv')
  # p535_total1 <- cowplot::plot_grid(res_SA535_UnRx$plg, p535_total, ncol=1, rel_heights = c(1,10), labels = c('SA535',''))
  png(paste0(output_dir,datatag,".png"), height = 2*500, width=2*1100,res = 2*72)
  print(p535_total)
  dev.off()
  
  res_SA535_RxH <- readRDS(paste0(input_dir,'figs/SA535_RxH_dlp_prevalence.rds'))
  res_SA535_Rx <- readRDS(paste0(input_dir,'figs/SA535_Rx_dlp_prevalence.rds'))
  res_SA535_UnRx <- readRDS(paste0(input_dir,'figs/SA535_UnRx_dlp_prevalence.rds'))
  res_dlp_SA535_lg <- readRDS(paste0(input_dir,'figs/',datatag,'_clone_legend_plt.rds'))
  res_dlp_SA535_lg$plg
  res_SA535_RxH$p
  res_SA535_Rx$p
  res_SA535_UnRx$p
  
} 
  
  # saveRDS(p535, paste0(output_dir, datatag, "_total_UMAPs.rds"))
  # return(p535)
# }

plot_SA609 <- function(output_dir){
  output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/figs_rna/'
  # dir.create(output_dir)
  input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/'
  
  datatag <- 'SA609'
  basename <- datatag
  clonealign_stat <- data.table::fread(paste0(dirname(output_dir),'/clonealign_labels.csv'))
  
  unique(clonealign_stat$unique_clone)
  clonealign_stat <- clonealign_stat %>%
    dplyr::filter(datatag==basename)%>%
    dplyr::select(cell_id, unique_clone)%>%
    dplyr::rename(clone=unique_clone)
  colnames(clonealign_stat)
  rm(umap_df)
  umap_df <- data.table::fread(paste0(input_dir,datatag,'_norm_umap.csv')) %>% as.data.frame()
  dim(umap_df)
  
  # umap_df$scell_id <- paste0(umap_df$id,'_',umap_df$Barcode)
  # sum(umap_df$cell_id %in% clonealign_stat$cell_id)
  # umap_df <- umap_df %>% left_join(clonealign_stat, by=c('cell_id'))
  
  cols_use <- get_color_clones(datatag)
  tps <- unique(umap_df$timepoint)
  rna_SA609 <- list()
  for(obs_treatment in c('UnRx','Rx','RxH')){
    # res <- viz_umap_obs_clones(umap_df, cols_use, datatag, 
    #                            output_dir, obs_treatment, NULL, 
    #                            obs_treatment, F)
    # rna_SA609[[paste0(obs_treatment,'_total')]] <- res
    for(obs_passage in tps){
      plottitle= paste0(obs_treatment,': ',obs_passage)
      res <- viz_umap_obs_clones(umap_df, cols_use, datatag, 
                                 output_dir, obs_treatment, obs_passage, 
                                 plottitle, F)
      rna_SA609[[paste0(obs_treatment,'_',obs_passage)]] <- res
    }
  }
  
  
  p609_total <- cowplot::plot_grid(rna_SA609$UnRx_X3$p, rna_SA609$UnRx_X4$p, rna_SA609$UnRx_X5$p, rna_SA609$UnRx_X6$p, rna_SA609$UnRx_X7$p,res_SA609_UnRx$p,
                                   res_dlp_SA609_lg$plg, rna_SA609$Rx_X4$p,rna_SA609$Rx_X5$p, rna_SA609$Rx_X6$p,rna_SA609$Rx_X7$p,res_SA609_Rx$p,
                                   NULL, NULL, rna_SA609$RxH_X5$p,rna_SA609$RxH_X6$p, rna_SA609$RxH_X7$p,res_SA609_RxH$p,
                                   ncol=6, align = 'hv')
  # p609_total1 <- cowplot::plot_grid(res_SA609_Rx$plg, p609_total, ncol=1, rel_heights = c(1,10), labels = c('SA609',''))
  png(paste0(output_dir,datatag,".png"), height = 2*500, width=2*800,res = 2*72)
  print(p609_total)
  dev.off()
  
  
  res_SA609_UnRx <- readRDS(paste0(input_dir,'figs/',datatag,'_UnRx_dlp_prevalence.rds'))
  res_SA609_Rx <- readRDS(paste0(input_dir,'figs/',datatag,'_Rx_dlp_prevalence.rds'))
  res_SA609_RxH <- readRDS(paste0(input_dir,'figs/',datatag,'_RxH_dlp_prevalence.rds'))
  res_dlp_SA609_lg <- readRDS(paste0(input_dir,'figs/',datatag,'_clone_legend_plt.rds'))
  res_dlp_SA609_lg$plg
}


plot_SA1035 <- function(output_dir){
  output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/figs_rna/'
  # dir.create(output_dir)
  input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/'
  
  datatag <- 'SA1035'
  basename <- datatag
  clonealign_stat <- data.table::fread(paste0(dirname(output_dir),'/clonealign_labels.csv'))
  
  # unique(clonealign_stat$unique_clone)
  clonealign_stat <- clonealign_stat %>%
    dplyr::filter(datatag==basename)%>%
    dplyr::select(cell_id, unique_clone)%>%
    dplyr::rename(clone=unique_clone)
  colnames(clonealign_stat)
  rm(umap_df)
  umap_df <- data.table::fread(paste0(input_dir,datatag,'_norm_umap.csv')) %>% as.data.frame()
  dim(umap_df)
  umap_df$cell_id[1]
  # umap_df$scell_id <- paste0(umap_df$id,'_',umap_df$Barcode)
  # sum(umap_df$cell_id %in% clonealign_stat$cell_id)
  # umap_df <- umap_df %>% left_join(clonealign_stat, by=c('cell_id'))
  
  cols_use <- get_color_clones(datatag)
  tps <- unique(umap_df$timepoint)
  print(tps)
  rna_SA1035 <- list()
  for(obs_treatment in c('UnRx','Rx','RxH')){
    # res <- viz_umap_obs_clones(umap_df, cols_use, datatag, 
    #                            output_dir, obs_treatment, NULL, 
    #                            obs_treatment, F)
    # rna_SA609[[paste0(obs_treatment,'_total')]] <- res
    for(obs_passage in tps){
      plottitle= paste0(obs_treatment,': ',obs_passage)
      res <- viz_umap_obs_clones(umap_df, cols_use, datatag, 
                                 output_dir, obs_treatment, obs_passage, 
                                 plottitle, F)
      rna_SA1035[[paste0(obs_treatment,'_',obs_passage)]] <- res
    }
  }
  
  # all SA1035 drug holiday are unassigned clones, double check input data
  p1035_total <- cowplot::plot_grid(rna_SA1035$UnRx_X4$p, rna_SA1035$UnRx_X5$p, 
                                    rna_SA1035$UnRx_X6$p, rna_SA1035$UnRx_X7$p, 
                                    rna_SA1035$UnRx_X8$p,res_SA1035_UnRx$p,  # untreated
                                    res_dlp_SA1035_lg$plg, rna_SA1035$Rx_X5$p, 
                                    rna_SA1035$Rx_X6$p,rna_SA1035$Rx_X7$p,
                                    rna_SA1035$Rx_X8$p,res_SA1035_Rx$p,   # treated
                                    NULL, NULL,
                                    rna_SA1035$RxH_X6$p, rna_SA1035$RxH_X7$p,
                                    rna_SA1035$RxH_X8$p,res_SA1035_RxH$p,
                                    ncol=6, align = 'hv')
  # p609_total1 <- cowplot::plot_grid(res_SA609_Rx$plg, p609_total, ncol=1, rel_heights = c(1,10), labels = c('SA609',''))
  png(paste0(output_dir,datatag,".png"), height = 2*500, width=2*800,res = 2*72)
  print(p1035_total)
  dev.off()
  
  
  res_SA1035_UnRx <- readRDS(paste0(input_dir,'figs/',datatag,'_UnRx_dlp_prevalence.rds'))
  res_SA1035_Rx <- readRDS(paste0(input_dir,'figs/',datatag,'_Rx_dlp_prevalence.rds'))
  res_SA1035_RxH <- readRDS(paste0(input_dir,'figs/',datatag,'_RxH_dlp_prevalence.rds'))
  res_dlp_SA1035_lg <- readRDS(paste0(input_dir,'figs/',datatag,'_clone_legend_plt.rds'))
  res_dlp_SA1035_lg$plg
  res_SA1035_UnRx$p
}


plot_addition_pdxs <- function(output_dir){
  
  output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/figs_rna/'
  # dir.create(output_dir)
  input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/'
  
  datatag <- 'SA501'
  basename <- datatag
  clonealign_stat <- data.table::fread(paste0(dirname(output_dir),'/clonealign_labels.csv'))
  
  # unique(clonealign_stat$unique_clone)
  clonealign_stat <- clonealign_stat %>%
    dplyr::filter(datatag==basename)%>%
    dplyr::select(cell_id, unique_clone)%>%
    dplyr::rename(clone=unique_clone)
  colnames(clonealign_stat)
  rm(umap_df)
  umap_df <- data.table::fread(paste0(input_dir,datatag,'_norm_umap.csv')) %>% as.data.frame()
  dim(umap_df)
  umap_df$cell_id[1]
  # umap_df$scell_id <- paste0(umap_df$id,'_',umap_df$Barcode)
  # sum(umap_df$cell_id %in% clonealign_stat$cell_id)
  # umap_df <- umap_df %>% left_join(clonealign_stat, by=c('cell_id'))
  
  cols_use <- get_color_clones(datatag)
  tps <- unique(umap_df$timepoint)
  print(tps)
  rna_SA501 <- list()
  obs_treatment <- 'UnRx'
  for(obs_passage in tps){
    plottitle= paste0(datatag, ' ', obs_treatment,': ',obs_passage)
    res <- viz_umap_obs_clones(umap_df, cols_use, datatag, 
                               output_dir, obs_treatment, obs_passage, 
                               plottitle, F)
    rna_SA501[[paste0(obs_treatment,'_',obs_passage)]] <- res
  }
  
  res_SA501_UnRx <- readRDS(paste0(input_dir,'figs/',datatag,'_dlp_prevalence.rds'))
  
  res_dlp_SA501_lg <- readRDS(paste0(input_dir,'figs/',datatag,'_clone_legend_plt.rds'))
  
  
  
  datatag <- 'SA530'
  basename <- datatag
  clonealign_stat <- data.table::fread(paste0(dirname(output_dir),'/clonealign_labels.csv'))
  
  # unique(clonealign_stat$unique_clone)
  clonealign_stat <- clonealign_stat %>%
    dplyr::filter(datatag==basename)%>%
    dplyr::select(cell_id, unique_clone)%>%
    dplyr::rename(clone=unique_clone)
  colnames(clonealign_stat)
  rm(umap_df)
  umap_df <- data.table::fread(paste0(input_dir,datatag,'_norm_umap.csv')) %>% as.data.frame()
  dim(umap_df)
  umap_df$cell_id[1]
  # umap_df$scell_id <- paste0(umap_df$id,'_',umap_df$Barcode)
  # sum(umap_df$cell_id %in% clonealign_stat$cell_id)
  # umap_df <- umap_df %>% left_join(clonealign_stat, by=c('cell_id'))
  
  cols_use <- get_color_clones(datatag)
  tps <- unique(umap_df$timepoint)
  print(tps)
  rna_SA530 <- list()
  obs_treatment <- 'UnRx'
  for(obs_passage in tps){
    plottitle= paste0(datatag, ' ', obs_treatment,': ',obs_passage)
    res <- viz_umap_obs_clones(umap_df, cols_use, datatag, 
                               output_dir, obs_treatment, obs_passage, 
                               plottitle, F)
    rna_SA530[[paste0(obs_treatment,'_',obs_passage)]] <- res
  }
  res_SA530_UnRx <- readRDS(paste0(input_dir,'figs/',datatag,'_dlp_prevalence.rds'))
  
  res_dlp_SA530_lg <- readRDS(paste0(input_dir,'figs/',datatag,'_clone_legend_plt.rds'))
  
  p501_530_total <- cowplot::plot_grid(rna_SA501$UnRx_X2$p, res_SA501_UnRx$p, res_dlp_SA501_lg$plg, 
                                    rna_SA530$UnRx_X3$p,res_SA530_UnRx$p,res_dlp_SA530_lg$plg,
                                    nrow=1)
  # p609_total1 <- cowplot::plot_grid(res_SA609_Rx$plg, p609_total, ncol=1, rel_heights = c(1,10), labels = c('SA609',''))
  png(paste0(output_dir,datatag,".png"), height = 2*180, width=2*800,res = 2*72)
  print(p501_530_total)
  dev.off()
  
  datatag <- 'SA604'
  basename <- datatag
  clonealign_stat <- data.table::fread(paste0(dirname(output_dir),'/clonealign_labels.csv'))
  
  # unique(clonealign_stat$unique_clone)
  clonealign_stat <- clonealign_stat %>%
    dplyr::filter(datatag==basename)%>%
    dplyr::select(cell_id, unique_clone)%>%
    dplyr::rename(clone=unique_clone)
  colnames(clonealign_stat)
  rm(umap_df)
  umap_df <- data.table::fread(paste0(input_dir,datatag,'_norm_umap.csv')) %>% as.data.frame()
  dim(umap_df)
  umap_df$cell_id[1]
  # umap_df$scell_id <- paste0(umap_df$id,'_',umap_df$Barcode)
  # sum(umap_df$cell_id %in% clonealign_stat$cell_id)
  clonealign_stat$cell_id[1]
  umap_df <- umap_df %>% left_join(clonealign_stat, by=c('cell_id'))
  umap_df$cell_id <- paste0(umap_df$id,'_', umap_df$Barcode)
  sum(umap_df$cell_id %in% clonealign_stat$cell_id)
  umap_df <- umap_df %>% left_join(clonealign_stat, by=c('cell_id'))
  
  cols_use <- get_color_clones(datatag)
  tps <- unique(umap_df$timepoint)
  print(tps)
  rna_SA604 <- list()
  obs_treatment <- 'UnRx'
  for(obs_passage in tps){
    plottitle= paste0(datatag, ' ', obs_treatment,': ',obs_passage)
    res <- viz_umap_obs_clones(umap_df, cols_use, datatag, 
                               output_dir, obs_treatment, obs_passage, 
                               plottitle, F)
    rna_SA604[[paste0(obs_treatment,'_',obs_passage)]] <- res
  }
  res_SA604_UnRx <- readRDS(paste0(input_dir,'figs/',datatag,'_UnRx_dlp_prevalence.rds'))
  res_SA604_UnRx$p
  res_dlp_SA604_lg <- readRDS(paste0(input_dir,'figs/',datatag,'_clone_legend_plt.rds'))
  res_dlp_SA604_lg$plg
  
  p604_total <- cowplot::plot_grid(rna_SA604$UnRx_X7$p, 
                                   rna_SA604$UnRx_X8$p, rna_SA604$UnRx_X9$p,
                                   res_SA604_UnRx$p,res_dlp_SA604_lg$plg,
                                   nrow=1)
  # p609_total1 <- cowplot::plot_grid(res_SA609_Rx$plg, p609_total, ncol=1, rel_heights = c(1,10), labels = c('SA609',''))
  png(paste0(output_dir,datatag,".png"), height = 2*180, width=2*800,res = 2*72)
  print(p604_total)
  dev.off()
  
} 

plot_total <- function(){
  saveRDS(p501_530_total, paste0(output_dir,"p501_530_total.rds"))
  saveRDS(p604_total, paste0(output_dir,"p604_total.rds"))
  saveRDS(p609_total, paste0(output_dir,"p609_total.rds"))
  saveRDS(p535_total, paste0(output_dir,"p535_total.rds"))
  saveRDS(p1035_total, paste0(output_dir,"p1035_total.rds"))
  
  p_all <- cowplot::plot_grid(p501_530_total,p604_total, 
                                   p609_total, p535_total, p1035_total,
                                   nrow=5, rel_heights = c(1,1,3,3,3))
  desc <- 'all_series_UMAPs'
  png(paste0(output_dir,desc,".png"), height = 2*110*11, width=2*900,res = 2*72)
  print(p_all)
  dev.off()
}

