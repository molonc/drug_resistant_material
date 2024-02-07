library(dplyr)
extract_copy_number <- function(){
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  # source(paste0(base_dir, 'scripts/drug_manuscript/cnv_viz_utils.R'))
  input_dir <- paste0(base_dir, 'materials/cell_clones/')
  datatag <- 'SA609'
  copynumber_fn <- paste0(input_dir, datatag, '_total_merged_filtered_states.csv.gz')
  copynumber <- data.table::fread(copynumber_fn) %>% as.data.frame()
  rownames(copynumber) <- copynumber$V1
  copynumber$V1 <- NULL
  dim(copynumber)
  rownames(copynumber)[1]
  
 
  df <- data.frame(cell_id=colnames(copynumber), sample_id=get_sample_id(colnames(copynumber))) 
  summary(as.factor(df$sample_id))
  obs_samples <- c('SA609X3XB01584','SA609X4XB03080','SA609X5XB03223',
                   'SA609X6XB03447','SA609X7XB03554','SA609X4XB003083',
                   'SA609X5XB03230','SA609X6XB03404','SA609X7XB03505',
                   'SA609X5XB03231','SA609X6XB03401','SA609X7XB03510')
  sum(df$sample_id %in% obs_samples)
  cell_clones <- data.table::fread(paste0(input_dir, 'SA609_cell_clones_total.csv.gz')) %>% as.data.frame()
  
  dim(cell_clones)  
  cell_clones$cell_id[1]
  summary(as.factor(cell_clones$sample))
  sum(cell_clones$sample %in% obs_samples)
  length(unique(cell_clones$sample))
  cell_clones$sample_id <- get_sample_id(cell_clones$cell_id)
  
  sids <- lapply(cell_clones$sample, function(x) {
    return(substr(x, nchar(x)-3, nchar(x)))
  })
  obs_samples_ids <- lapply(obs_samples, function(x) {
    return(substr(x, nchar(x)-3, nchar(x)))
  })
  sids <- unlist(sids)
  obs_samples_ids <- unlist(obs_samples_ids)
  sum(sids %in% obs_samples_ids)
  
  cell_clones <- cell_clones %>%
    dplyr::filter(sample_id %in% obs_samples)
  dim(cell_clones)  
  summary(as.factor(cell_clones$clone_id))
  summary(as.factor(cell_clones$clone_id))
  cell_clones <- cell_clones %>%
    dplyr::filter(clone_id=='R')
  unique(cell_clones$sample)
}  
##---------------------------------------------------------------------------------------
base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
source(paste0(base_dir, 'scripts/drug_manuscript/viz_umap_figs/viz_umaps.R'))
source(paste0(base_dir, 'scripts/drug_manuscript/cnv_viz_utils.R'))

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
  datatag <- 'SA609'
  series <- c('SA609')
  pts_lb <- c('Pt4') 
  names(pts_lb) <- series
  
  dlp_heatmap_ls <- list()
  for(datatag in series){
    # copynumber_fn <- paste0(input_dir, datatag, '_total_merged_filtered_states.csv.gz')
    copynumber_fn <- paste0(input_dir, datatag, '_total_merged_filtered_states_longformat.csv.gz')
    cellclone_fn <- paste0(input_dir, datatag, '_cell_clones.csv.gz')
    library_grouping_fn <- paste0(input_dir, datatag, '_DLP_library_groupings.csv')
    df_cnv <- get_median_genotype_v3(copynumber_fn, datatag, save_dir,
                                     cellclone_fn, library_grouping_fn) 
    print(dim(df_cnv))
    df_cnv <- df_cnv %>%
      dplyr::rename(cnv=copy_number)
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
  coef_df <- data.table::fread(paste0(base_dir,'materials/fitness_paper_DLP/sa609_extended_data_fig8.csv.gz')) %>% as.data.frame()#s_coeff_filter.csv.gz
  coef_df$datatag <- 'SA609'
  coef_df$colour <- NULL
  ######### Untreated ######### 
  rm(dat)
  coef_df$s <- coef_df$s - 1
  dat <- coef_df %>%
    dplyr::filter(label=="Left (-Rx Line 2: SA609)")                  
  dim(dat)
  
  tag <- 'SA609'
  plottitle <- 'UnRx'
  p_coeff_UnRx <- plot_box_posterior_v2(dat, color_fn, tag, plottitle, ymax=NULL, ymin=NULL,   
                       clone_axis_dir = 'bottom', sort_by = 'median', add_linebreak = FALSE)
    
  p_coeff_UnRx
  saveRDS(p_coeff_UnRx, paste0(output_dir, tag, "_UnRx_coeff_fitness.rds"))
  
  ######### Treated ######### 
  rm(dat)
  dat <- coef_df %>%
    dplyr::filter(label=="Right (+Rx Line 1: SA609)")                  
  dim(dat)
  
  tag <- 'SA609'
  plottitle <- 'Rx'
  p_coeff_Rx <- plot_box_posterior_v2(dat, color_fn, tag, plottitle, ymax=1.4, ymin=0.6,   
                             clone_axis_dir = 'bottom', sort_by = 'median', add_linebreak = FALSE)
  p_coeff_Rx  
  saveRDS(p_coeff_Rx, paste0(output_dir, tag, "_Rx_coeff_fitness.rds"))
  
}

plot_10x_UMAP <- function()
{  
  
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  source(paste0(base_dir, 'scripts/drug_manuscript/viz_umap_figs/viz_umaps.R'))
  
  input_dir <- paste0(base_dir, 'materials/umap_figs/')
  output_dir <- paste0(base_dir, 'materials/umap_figs/figs_rna/')
  
  datatag <- 'SA609'
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
  dim(umap_df)
  
  # umap_df$scell_id <- paste0(umap_df$id,'_',umap_df$Barcode)
  # sum(umap_df$cell_id %in% clonealign_stat$cell_id)
  # umap_df <- umap_df %>% left_join(clonealign_stat, by=c('cell_id'))
  
  # color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v3.csv.gz')
  color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v4.csv.gz')
  cols_use <- get_color_clones(datatag, color_fn) # predefined clone color code for DLP and 10x inferred clones clonealign
  
  unique(umap_df$clone)
  umap_df <- umap_df %>%
    dplyr::mutate(clone = get_unique_clone_id(clone)) %>%
    dplyr::mutate(clone = case_when(clone=='R' ~ 'A',
                                       TRUE  ~ clone))
  dim(umap_df)
  
  ## For pseudotime clonal labels, revision manuscript
  # data.table::fwrite(umap_df, '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/clone_labels_unique_SA609.csv.gz')
  
  
  res_prop10x_SA609 <- plot_fill_barplot_wholedataset_rnaseq(umap_df, cols_use, output_dir, 
                                                             datatag, plottitle=NULL, plotlegend=F, facet_order='vertical')
  
  # dev.off()
  # res_prop10x_SA609$p
  saveRDS(res_prop10x_SA609, paste0(output_dir, datatag, "_prevalence_clone_clonealign10x.rds"))
  
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
  saveRDS(rna_SA609, paste0(output_dir,datatag,"_umap_plt.rds"))
  # res_SA609_UnRx$p,res_dlp_SA609_lg$plg, 
  # res_SA609_Rx$p,res_SA609_RxH$p,
  # rna_SA609$UnRx_X3$p, 
  # clone_plg_609 <- plot_clone_color_legend('SA609', base_dir, ncols_grid=2)
  
  # p609_total <- cowplot::plot_grid(rna_SA609$UnRx_X4$p, rna_SA609$UnRx_X5$p, rna_SA609$UnRx_X6$p, rna_SA609$UnRx_X7$p,
  #                                  rna_SA609$Rx_X4$p,rna_SA609$Rx_X5$p, rna_SA609$Rx_X6$p,rna_SA609$Rx_X7$p,
  #                                  clone_plg_609, rna_SA609$RxH_X5$p,rna_SA609$RxH_X6$p, rna_SA609$RxH_X7$p,
  #                                  ncol=4, align = 'hv')
  
 
  
   
}


plot_DLP_barplot <- function(){
  
  ##---------------------------------------------------------------------------------------
  ##-----------------------------------SA604 barplot prevalence from sitka tree -----------
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  input_dir <- paste0(base_dir, 'materials/cell_clones/')
  output_dir <- paste0(base_dir, 'materials/umap_figs/figs_rna/')
  # source(paste0(base_dir, 'scripts/drug_manuscript/viz_umap_figs/viz_umaps.R'))  
  datatag <- 'SA609'
  library_grouping_fn <- paste0(input_dir, datatag, '_DLP_library_groupings.csv.gz')
  cell_clones_fn <- paste0(input_dir, datatag, '_cell_clones_total.csv.gz')
  cell_clones <- get_meta_data(cell_clones_fn, library_grouping_fn, datatag)
  unique(cell_clones$treatment_desc)
  table(cell_clones$treatmentSt, cell_clones$timepoint)
  
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
  res_barDLP_609 <- plot_fill_barplot_wholedataset_v2(cell_clones, cols_use, output_dir, 
                                                      datatag, plottitle=NULL, plotlegend=F, facet_order='vertical')
  # res_barDLP_609$p
  # res_barDLP_609$plg
  
  saveRDS(res_barDLP_609, paste0(output_dir, datatag, "_barplot_DLP.rds"))
  
  
}



## Eric
plot_SUPP_fig2 <- function(){
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  input_dir <- paste0(base_dir, 'materials/cell_clones/')
  output_dir <- paste0(base_dir, 'materials/umap_figs/figs_rna/')
  source(paste0(base_dir, 'scripts/drug_manuscript/viz_umap_figs/viz_umaps.R'))
  
  datatag <- 'SA609'
  # medianCNV_plt_topplot <- res$p
  medianCNV_plt_topplot <- readRDS(paste0(output_dir, datatag, '_DLP_heatmap.rds'))
  medianCNV_plt_topplot1 <- cowplot::plot_grid(medianCNV_plt_topplot$cnv_plot, 
                                               medianCNV_plt_topplot$plg, rel_heights = c(1,0.1), ncol=1,
                                               labels=c('  a ',NULL),label_size=11, hjust=0)
  
  # quantiles of selection coefficients of fitClone model 
  # fits to each branch with respect to the reference Clone C in TNBC-SA609
  p_coeff_UnRx <- readRDS(paste0(output_dir, datatag, "_UnRx_coeff_fitness.rds"))
  p_coeff_Rx <- readRDS(paste0(output_dir, datatag, "_Rx_coeff_fitness.rds"))
  plt_topplot2 <-cowplot::plot_grid(NULL, p_coeff_UnRx, p_coeff_Rx, ncol=1, 
                                    rel_heights = c(0.2, 1, 1), labels = c('  b Fitness coefficients',NULL,NULL),
                                    label_size=11, hjust=0)
  
  p609_top_plt <- cowplot::plot_grid(medianCNV_plt_topplot1, plt_topplot2, rel_widths = c(3.8,2.2), nrow=1)#+ 
    # theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  clone_plg_609 <- plot_clone_color_legend('SA609', base_dir, ncols_grid=2)
  clone_plg_609
  
  
  res_barDLP_609 <- readRDS(paste0(output_dir, datatag, "_barplot_DLP.rds"))
  # res_barDLP_609$p
  res_prop10x_SA609 <- readRDS(paste0(output_dir, datatag, "_prevalence_clone_clonealign10x.rds"))
  rna_SA609 <- readRDS(paste0(output_dir,datatag,"_umap_plt.rds"))
  
  p609_leftside <- cowplot::plot_grid(rna_SA609$UnRx_X4$p, rna_SA609$UnRx_X5$p, rna_SA609$UnRx_X6$p, rna_SA609$UnRx_X7$p,
                                   rna_SA609$Rx_X4$p,rna_SA609$Rx_X5$p, rna_SA609$Rx_X6$p,rna_SA609$Rx_X7$p,
                                   clone_plg_609, rna_SA609$RxH_X5$p,rna_SA609$RxH_X6$p, rna_SA609$RxH_X7$p,
                                   ncol=4, align = 'hv')#+
    # theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  p_prevalence <- cowplot::plot_grid(res_prop10x_SA609$p, res_barDLP_609$p, rel_widths = c(1, 1), nrow=1)
  p_prevalence1 <- cowplot::plot_grid(NULL, p_prevalence, rel_heights = c(0.05, 1), ncol=1)
  p609_bottom <- cowplot::plot_grid(p609_leftside, NULL, p_prevalence1,
                                    rel_widths = c(4, 0.05, 2), nrow=1)#+ 
    # theme(plot.background = element_rect(fill = "white", colour = "white"))
  p609_bottom1 <- cowplot::plot_grid(NULL, NULL, NULL,labels=c('  c 10x UMAP across treatments & time points',
                                                               '  d 10x inferred \nclones prevalence','   e DLP clones\n prevalence'),
                                     label_size=11, hjust=0, rel_widths = c(4,1,1),nrow=1)
  p609_total <- cowplot::plot_grid(p609_top_plt, p609_bottom1, p609_bottom, rel_heights = c(2,0.2,3), ncol=1) + 
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  # ggsave(paste0(output_dir,"SUPPFig2_SA609.pdf"),
  #        plot = p609_total,
  #        height = 11,
  #        width = 11,
  #        useDingbats=F, # from Sam's suggestion
  #        dpi = 150
  # )
  
  ggsave(paste0(output_dir,"SUPPFig2_SA609.svg"),
         plot = p609_total,
         height = 12,
         width = 11,
         # type = "cairo-png",
         dpi=200
  )
  
}








# dat: data.frame with columns: (id, s)
# clone_dic: a data.frame with columns: (id, K, letters, pretty_names, colour)
# plot_box_posterior <- function(dat, clone_dic, show_title=TRUE, ref_clone_letter = '', clone_axis_dir = 'bottom', sort_by = 'median', add_linebreak = FALSE) {
#   
#   clone_dic$K <- as.character(clone_dic$K)
#   dat <- dplyr::right_join(dat, clone_dic)
#   dat <- dat[!is.na(dat$id), ]
#   
#   
#   # Sort clones by their posterior mean fitness values
#   if (sort_by == 'median') {
#     s_means <- dat %>% dplyr::group_by(letters) %>% dplyr::summarise(mean_s = median(s)) %>% dplyr::arrange(mean_s)
#   } else {
#     s_means <- dat %>% dplyr::group_by(letters) %>% dplyr::summarise(mean_s = mean(s)) %>% dplyr::arrange(mean_s)
#   }
#   
#   dat$letters <- factor(x = dat$letters, levels = s_means$letters, ordered = T)
#   
#   
#   myColors <- clone_dic$colour
#   names(myColors) <- clone_dic$letters
#   
#   dat$pretty_names <- factor(dat$pretty_names, levels = names(myColors), ordered = T)
#   stopifnot(all(levels(dat$pretty_names) == as.character(names(myColors))))
#   
#   # change to s + 1
#   dat$s <- dat$s + 1
#   
#   # Set the relevant s-values, length of the whiskers  
#   # Put the min and max in range of the whiskers...
#   whisker_coef <- 1.5
#   
#   lims <- dat %>% dplyr::group_by(K) %>% dplyr::summarise(iqr = IQR(s), median = median(s), q1 = quantile(s, .25), q3 = quantile(s, .75), bottom = min(s), top = max(s)) 
#   
#   p <- 
#     dat %>% 
#     ggplot(aes(letters, s, color = letters, fill = letters)) +
#     geom_boxplot(outlier.shape = NA, coef = whisker_coef) + 
#     stat_boxplot(geom = 'errorbar', width = 0.3, size = .5, coef = whisker_coef) + 
#     scale_x_discrete(limits = levels(dat$letters), breaks = levels(dat$letters), labels = levels(dat$letters), position = clone_axis_dir) + 
#     coord_flip(ylim = c(min(lims$bottom), max(lims$top))) + 
#     geom_hline(yintercept = 1.0, linetype = 'longdash', size = .5) +
#     scale_color_manual(values = myColors) + 
#     scale_fill_manual(values = myColors) + 
#     xlab('') + ylab(sprintf('1 + s (relative to clone %s)', ref_clone_letter)) +
#     theme_light(base_size = 20) + 
#     fitclone_get_theme_no_grid() + 
#     theme(axis.ticks.y = element_blank())
#   
#   if (add_linebreak) {
#     p <- p + ylab(sprintf('1 + s\n(rel. to clone %s)', ref_clone_letter))
#   }
#   
#   
#   if (show_title) {
#     p <- p + ggtitle(paste0(datatag), subtitle = get_ID_from_edge_list_path(edge_list_path))
#   } else {
#     p <- p + theme(legend.position = "none")
#   }
#   
#   data <- ggplot_build(p)$data[[1]]
#   p <- p + geom_segment(inherit.aes = FALSE, data=data, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", size=.5)
#   
#   return(p)
# }


# plot_box_posterior_v2 <- function(dat, clone_dic, show_title=TRUE, ref_clone_letter = '', 
#                                   clone_axis_dir = 'bottom', sort_by = 'median', add_linebreak = FALSE) {
#   
#   letters=clone_name
#   s=median_1_plus_s
#   pretty_names
#   coeff_df$clone_name
#   
#   # Sort clones by their posterior mean fitness values
#   if (sort_by == 'median') {
#     s_means <- dat %>% dplyr::group_by(letters) %>% dplyr::summarise(mean_s = median(s)) %>% dplyr::arrange(mean_s)
#   } else {
#     s_means <- dat %>% dplyr::group_by(letters) %>% dplyr::summarise(mean_s = mean(s)) %>% dplyr::arrange(mean_s)
#   }
#   
#   dat$letters <- factor(x = dat$letters, levels = s_means$letters, ordered = T)
#   
#   
#   myColors <- clone_dic$colour
#   names(myColors) <- clone_dic$letters
#   
#   dat$pretty_names <- factor(dat$pretty_names, levels = names(myColors), ordered = T)
#   stopifnot(all(levels(dat$pretty_names) == as.character(names(myColors))))
#   
#   # change to s + 1
#   # dat$s <- dat$s + 1
#   
#   # Set the relevant s-values, length of the whiskers  
#   # Put the min and max in range of the whiskers...
#   whisker_coef <- 1.5
#   
#   lims <- dat %>% dplyr::group_by(K) %>% dplyr::summarise(iqr = IQR(s), median = median(s), q1 = quantile(s, .25), q3 = quantile(s, .75), bottom = min(s), top = max(s)) 
#   
#   p <- 
#     dat %>% 
#     ggplot(aes(letters, s, color = letters, fill = letters)) +
#     geom_boxplot(outlier.shape = NA, coef = whisker_coef) + 
#     stat_boxplot(geom = 'errorbar', width = 0.3, size = .5, coef = whisker_coef) + 
#     scale_x_discrete(limits = levels(dat$letters), breaks = levels(dat$letters), labels = levels(dat$letters), position = clone_axis_dir) + 
#     coord_flip(ylim = c(min(lims$bottom), max(lims$top))) + 
#     geom_hline(yintercept = 1.0, linetype = 'longdash', size = .5) +
#     scale_color_manual(values = myColors) + 
#     scale_fill_manual(values = myColors) + 
#     xlab('') + ylab(sprintf('1 + s (relative to clone %s)', ref_clone_letter)) +
#     theme_light(base_size = 20) + 
#     fitclone_get_theme_no_grid() + 
#     theme(axis.ticks.y = element_blank())
#   
#   if (add_linebreak) {
#     p <- p + ylab(sprintf('1 + s\n(rel. to clone %s)', ref_clone_letter))
#   }
#   
#   
#   p <- p + theme(legend.position = "none")
#   
#   data <- ggplot_build(p)$data[[1]]
#   p <- p + geom_segment(inherit.aes = FALSE, data=data, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", size=.5)
#   
#   return(p)
# }

# plot_coefficients <- function(){
#   base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
#   input_dir <- paste0(base_dir, 'materials/fitness_paper_DLP/')
#   output_dir <- paste0(base_dir, 'materials/umap_figs/figs_rna/')
#   coeff_df <- data.table::fread(paste0(input_dir, 'SUPP_Table2_fitness_coefficients.csv')) %>% as.data.frame()
#   
#   
#   coeff_df <- coeff_df %>%
#     dplyr::filter(datasetname=='SA609')
#   
# }

