
generate_DLP_collapsed_trees <- function(){
  script.basename <- '/home/htran/Projects/farhia_project/drug_resistant_material/scripts/corrupt_tree/src/tree_viz'
  source(paste0(script.basename, "/utils.R"))
  # source(paste0(script.basename, "/trajectory_utils_v2.R"))
  source(paste0(script.basename, "/ebola_tree_utils.R"))
  source(paste0(script.basename, "/collapse_tree_utils.R"))
  
  datatag <- 'SA501'
  results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA501_v2/'
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/'
  colorcode_fn <- paste0(base_dir,'umap_figs/colorcode_total_v4.csv')
  output_fn <- paste0(results_dir, 'tree_viz_dream/collapse_tree_wholedata.png')
  cellclones_fn <- paste0(results_dir, 'cell_clones.csv')
  newick_fn <- paste0(results_dir, 'tree.newick')
  # Note: clone id: R--> H
  collapse_tree_501 <- plot_ebola_tree_condition_main_fig2(datatag, results_dir, colorcode_fn,
                                                    output_fn, cellclones_fn, newick_fn)
  
  
  # collapse_tree_501
  
  datatag <- 'SA530'
  results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA530_v2/'
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/'
  colorcode_fn <- paste0(base_dir,'umap_figs/colorcode_total_v4.csv')
  output_fn <- paste0(results_dir, 'tree_viz_dream/collapse_tree_wholedata.png')
  cellclones_fn <- paste0(results_dir, 'cell_clones.csv')
  newick_fn <- paste0(results_dir, 'tree.newick')
  collapse_tree_530 <- plot_ebola_tree_condition_main_fig2(datatag, results_dir, colorcode_fn,
                                                    output_fn, cellclones_fn, newick_fn)
  
  
  # collapse_tree_530
  
  
  datatag <- 'SA604'
  results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA604_v2/'
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/'
  colorcode_fn <- paste0(base_dir,'umap_figs/colorcode_total_v4.csv')
  output_fn <- paste0(results_dir, 'tree_viz_dream/collapse_tree_wholedata.png')
  # Note: clone id: R1--> K, R2 --> K
  cellclones_fn <- paste0(results_dir, 'cell_clones.csv')
  newick_fn <- paste0(results_dir, 'tree.newick')
  collapse_tree_604 <- plot_ebola_tree_condition_main_fig2(datatag, results_dir, colorcode_fn,
                                                    output_fn, cellclones_fn, newick_fn)
  
  
  # collapse_tree_604
  
  
  datatag <- 'SA609'
  results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA609/'
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/'
  colorcode_fn <- paste0(base_dir,'umap_figs/colorcode_total_v4.csv')
  output_fn <- paste0(results_dir, 'tree_viz_dream/collapse_tree_wholedata.png')
  cellclones_fn <- paste0(results_dir, 'cell_clones.csv')
  newick_fn <- paste0(results_dir, 'SA609_omnibus.newick')
  collapse_tree_609 <- plot_ebola_tree_condition_main_fig2(datatag, results_dir, colorcode_fn,
                                                    output_fn, cellclones_fn, newick_fn)
  
  
  # collapse_tree_609
  
  datatag <- 'SA535'
  results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_fitness/'
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/'
  colorcode_fn <- paste0(base_dir,'umap_figs/colorcode_total_v4.csv')
  output_fn <- paste0(results_dir, 'tree_viz_dream/collapse_tree_wholedata.png')
  cellclones_fn <- paste0(results_dir, 'cell_clones.csv')
  newick_fn <- paste0(results_dir, 'SA535.newick')
  collapse_tree_535 <- plot_ebola_tree_condition_main_fig2(datatag, results_dir, colorcode_fn,
                                                           output_fn, cellclones_fn, newick_fn)
  
  
  # collapse_tree_535
  
  datatag <- 'SA1035'
  results_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_Tyler_v2/'
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/'
  colorcode_fn <- paste0(base_dir,'umap_figs/colorcode_total_v4.csv')
  output_fn <- paste0(results_dir, 'tree_viz_dream/collapse_tree_wholedata.png')
  cellclones_fn <- paste0(results_dir, 'cell_clones.csv')
  newick_fn <- paste0(results_dir, 'tree.newick')
  collapse_tree_1035 <- plot_ebola_tree_condition_main_fig2(datatag, results_dir, colorcode_fn,
                                                           output_fn, cellclones_fn, newick_fn)
  
  
  # collapse_tree_1035
  
  
}

base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/'
source(paste0(dirname(base_dir),'/scripts/drug_manuscript/viz_umap_figs/viz_umaps.R'))
source(paste0(dirname(base_dir),'/scripts/drug_manuscript/cnv_viz_utils.R'))
# base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/'

output_dir <- paste0(base_dir,'umap_figs/main_fig2/')

# Clone color definition
# color_df <- data.table::fread(paste0(base_dir,'umap_figs/colorcode_total_v3.csv'))
color_df <- data.table::fread(paste0(base_dir,'umap_figs/colorcode_total_v4.csv'))
# dim(color_df)
# summary(as.factor(color_df$datatag))

# clone_id, timepoint, treatmentstr
## Testing
#Noted: Remove X3 SA609, X5 SA535, X4 SA1035 early passage from analysis
plot_DLP_barplot_SA609 <- function(){
  
  ##---------------------------------------------------------------------------------------
  ##-----------------------------------SA604 barplot prevalence from sitka tree -----------
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  input_dir <- paste0(base_dir, 'materials/cell_clones/')
  output_dir <- paste0(base_dir, 'materials/umap_figs/main_fig2/')
  datatag <- 'SA609'
  library_grouping_fn <- paste0(input_dir, datatag, '_DLP_library_groupings.csv.gz')
  cell_clones_fn <- paste0(input_dir, datatag, '_cell_clones_total.csv.gz')
  cell_clones <- get_meta_data(cell_clones_fn, library_grouping_fn, datatag)
  unique(cell_clones$treatment_desc)
  table(cell_clones$treatmentSt, cell_clones$timepoint)
  
  summary(as.factor(cell_clones$clone_id))
  # cell_clones <- cell_clones %>%
  #   dplyr::filter(timepoint!='X3')
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
  color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v4.csv.gz')
  cols_use <- get_color_clones(datatag, color_fn) # predefined clone color code for DLP and 10x inferred clones clonealign
  
  
  cell_clones$treatment_desc <- paste0('Pt4-',cell_clones$treatment_desc)
  res_barDLP_609 <- plot_fill_barplot_wholedataset_v2(cell_clones, cols_use, output_dir, 
                                                      datatag, plottitle=NULL, plotlegend=F, 
                                                      facet_order='vertical', plot_yaxis=TRUE, yaxis_title="DLP Clone Prevalence")
  # res_barDLP_609$p
  # res_barDLP_609$plg
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  saveRDS(res_barDLP_609, paste0(output_dir, datatag, "_barplot_DLP_v2.rds"))
  
  return(res_barDLP_609)
}
res_barDLP_609 <- plot_DLP_barplot_SA609()


plot_DLP_barplot_SA535 <- function(){
  ##---------------------------------------------------------------------------------------
  ##-----------------------------------SA604 barplot prevalence from sitka tree -----------
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  input_dir <- paste0(base_dir, 'materials/cell_clones/')
  output_dir <- paste0(base_dir, 'materials/umap_figs/main_fig2/')
  
  datatag <- 'SA535'
  library_grouping_fn <- paste0(input_dir, datatag, '_DLP_library_groupings.csv.gz')
  cell_clones_fn <- paste0(input_dir, datatag, '_cell_clones.csv.gz')
  cell_clones <- get_meta_data(cell_clones_fn, library_grouping_fn, datatag)
  table(cell_clones$treatmentSt, cell_clones$timepoint)
  # cell_clones <- cell_clones %>%
  #   dplyr::filter(timepoint!='X5')
  # summary(as.factor(cell_clones$clone_id))
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
  color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v4.csv.gz')
  cols_use <- get_color_clones(datatag, color_fn) # predefined clone color code for DLP and 10x inferred clones clonealign
  
  
  cell_clones$treatment_desc <- paste0('Pt5-',cell_clones$treatment_desc)
  res_barDLP_535 <- plot_fill_barplot_wholedataset_v2(cell_clones, cols_use, output_dir, 
                                                      datatag, plottitle=NULL, plotlegend=F, facet_order='vertical')
  # res_barDLP_535$p
  # res_barDLP_535$plg
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  saveRDS(res_barDLP_535, paste0(output_dir, datatag, "_barplot_DLP_v2.rds"))
  return(res_barDLP_535)
  
}

res_barDLP_535 <- plot_DLP_barplot_SA535()

plot_DLP_barplot_SA1035 <- function(){
  ##---------------------------------------------------------------------------------------
  ##-----------------------------------SA604 barplot prevalence from sitka tree -----------
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  input_dir <- paste0(base_dir, 'materials/cell_clones/')
  output_dir <- paste0(base_dir, 'materials/umap_figs/main_fig2/')
  
  datatag <- 'SA1035'
  library_grouping_fn <- paste0(input_dir, datatag, '_DLP_library_groupings.csv.gz')
  cell_clones_fn <- paste0(input_dir, datatag, '_cell_clones.csv.gz')
  cell_clones <- get_meta_data(cell_clones_fn, library_grouping_fn, datatag)
  unique(cell_clones$treatment_desc)
  
  summary(as.factor(cell_clones$clone_id))
  table(cell_clones$treatmentSt, cell_clones$timepoint)
  # cell_clones <- cell_clones %>%
  #   dplyr::filter(timepoint!='X4')
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
  color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v4.csv.gz')
  cols_use <- get_color_clones(datatag, color_fn) # predefined clone color code for DLP and 10x inferred clones clonealign
  
  
  cell_clones$treatment_desc <- paste0('Pt6-',cell_clones$treatment_desc)
  res_barDLP_1035 <- plot_fill_barplot_wholedataset_v2(cell_clones, cols_use, output_dir, 
                                                       datatag, plottitle=NULL, plotlegend=F, facet_order='vertical')
  # res_barDLP_1035$p
  # res_barDLP_1035$plg
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  saveRDS(res_barDLP_1035, paste0(output_dir, datatag, "_barplot_DLP_v2.rds"))
  return(res_barDLP_1035)
  
}

res_barDLP_1035 <- plot_DLP_barplot_SA1035()

plot_tumour_growth_status <- function(output_dir){
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/'
  output_dir <- paste0(base_dir,'umap_figs/main_fig2/')
  input_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/tumour_growth_info/'
  status <- c('Progressive Disease','Stable Disease','Partial Response',
              'Complete Response','undefined')
  color_codes <- c('#E97451','#FAD5A5','lightgreen','#00FFFF','white')
  names(color_codes) <- status
  my_font <- "Helvetica"
  # library(dplyr)
  data <- data.table::fread(paste0(input_dir,'treatment_response_tumour_growth_v2.csv.gz'))

  plt_ls <- list()
  for(datatag in unique(data$Patient)){
    df <- data %>%
      dplyr::filter(Patient==datatag) %>%
      dplyr::select(Timepoint, Treatment, Response)
    # colnames(data)
    # View(data)
    # unique(df$Timepoint)
    if(datatag=='Pt4'){
      early_passage <- 'X3'
    }else if(datatag=='Pt5'){
      early_passage <- 'X5'
    }else{
      early_passage <- 'X4'
    }
    
    early_time <- data.frame(Timepoint=rep(early_passage,3),
                             Treatment=c('UnRx','RxH','Rx'),
                             Response='undefined')
    df <- dplyr::bind_rows(df, early_time)
    df$Treatment <- factor(df$Treatment, levels=c('RxH','Rx','UnRx'))
    df$Timepoint <- factor(df$Timepoint, levels=gtools::mixedsort(unique(df$Timepoint)))
    
    p <- ggplot(df, aes(x=Timepoint, y=Treatment, fill= Response, alpha=0.7)) + 
      geom_tile(colour = "white", size=1.2) + 
      # facet_wrap(Y~., scales = "free") + 
      theme_bw() + 
      scale_fill_manual(values = color_codes) + 
      theme(axis.line=element_blank(),
            axis.title = element_blank(), 
            # axis.title.x = element_blank(), 
            # axis.title.y = element_text(size=11), 
            # axis.text.y = element_text(size=11),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(size=12, family = my_font, colour = 'black'),
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            panel.grid = element_blank(),
            panel.border = element_blank())+ 
      labs(y=' ') #'- treatment status'paste0(datatag, '')
    plt_ls[[datatag]] <- p + theme(legend.position = 'none')
    
    saveRDS(p, paste0(output_dir, datatag, '_tumour_growth_status.rds'))
    
  }
  p1 <- p + guides(fill = guide_legend(title='Disease status',nrow = 1, size=4), alpha=NULL) + 
    theme(legend.position = 'bottom')
  plt_ls[['lg']] <- cowplot::ggdraw() + cowplot::draw_plot(cowplot::get_legend(p1))
  saveRDS(plt_ls[['lg']], paste0(output_dir, 'tumour_growth_legend.rds'))
  # total_tg <- tibble::tibble()
  # for(datatag in unique(data$Patient)){
  #   df <- data %>%
  #     dplyr::filter(Patient==datatag) %>%
  #     dplyr::select(Timepoint, Treatment, Response)
  #   # colnames(data)
  #   # View(data)
  #   # unique(df$Timepoint)
  #   if(datatag=='Pt4'){
  #     early_passage <- 'X3'
  #   }else if(datatag=='Pt5'){
  #     early_passage <- 'X5'
  #   }else{
  #     early_passage <- 'X4'
  #   }
  #   
  #   early_time <- data.frame(Timepoint=rep(early_passage,3),
  #                            Treatment=c('UnRx','RxH','Rx'),
  #                            Response='undefined')
  #   df <- dplyr::bind_rows(df, early_time)
  #   df$datatag <- datatag
  #   total_tg <- dplyr::bind_rows(total_tg, df)
  # }  
  # df <- total_tg
  # df$Treatment <- factor(df$Treatment, levels=c('RxH','Rx','UnRx'))
  # df$Timepoint <- factor(df$Timepoint, levels=gtools::mixedsort(unique(df$Timepoint)))
  # 
  # p <- ggplot(df, aes(x=Timepoint, y=Treatment, fill= Response, alpha=0.7)) + 
  #   geom_tile(colour = "white", size=1.2) + 
  #   facet_wrap(. ~ datatag, scales = "free") +
  #   theme_bw() + 
  #   scale_fill_manual(values = color_codes) + 
  #   theme(axis.line=element_blank(),
  #         axis.title = element_blank(), 
  #         # axis.title.x = element_blank(), 
  #         # axis.title.y = element_text(size=11), 
  #         # axis.text.y = element_text(size=11),
  #         axis.ticks.y = element_blank(),
  #         axis.text.y = element_blank(),
  #         axis.text.x = element_text(size=12, family = my_font),
  #         panel.background = element_rect(fill = "transparent", colour = NA),
  #         plot.background = element_rect(fill = "transparent", colour = NA),
  #         panel.grid = element_blank(),
  #         panel.border = element_blank(),
  #         legend.position = 'none'
  #   )
  # p
}
plot_mainFig2 <- function(){
  ## Loading DLP tree
  ## Noted: copy files into folder: "/home/htran/Projects/farhia_project/drug_resistant_material/materials/umap_figs/main_fig2/"
  # input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA501_v2/tree_viz_dream/'
  # dlp_501 <- readRDS(paste0(input_dir,'summary_tree.rds'))
  # 
  # input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA530_v2/tree_viz_dream/'
  # dlp_530 <- readRDS(paste0(input_dir,'summary_tree.rds'))
  # 
  # input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA604_v2/tree_viz_dream/'
  # dlp_604 <- readRDS(paste0(input_dir,'summary_tree.rds'))
  # 
  # input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA1035_new_encode/SA1035_Tyler_v2/tree_viz_dream/'
  # dlp_1035 <- readRDS(paste0(input_dir,'summary_tree.rds'))
  # 
  # input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_fitness/tree_viz_dream/'
  # dlp_535 <- readRDS(paste0(input_dir,'summary_tree.rds'))
  # 
  # input_dir <- '/home/htran/storage/datasets/drug_resistance/dlp_results/SA609/tree_viz_dream/'
  # dlp_609 <- readRDS(paste0(input_dir,'summary_tree.rds'))
  
  ## Loading DLP tree
  output_dir <- "/home/htran/Projects/farhia_project/drug_resistant_material/materials/umap_figs/main_fig2/"
  series <- c('SA501','SA530','SA604','SA609','SA535','SA1035')
  collapse_trees <- list()
  for(datatag in series){
    collapse_trees[[datatag]] <- readRDS(paste0(output_dir, datatag,'_summary_tree.rds'))
  }
  res_barDLP_609 <- readRDS(paste0(output_dir, "SA609_barplot_DLP_v2.rds"))
  res_barDLP_535 <- readRDS(paste0(output_dir, "SA535_barplot_DLP_v2.rds"))
  res_barDLP_1035 <- readRDS(paste0(output_dir, "SA1035_barplot_DLP_v2.rds"))
  
  ## Loading tumour growth status 
  tg_609 <- readRDS(paste0(output_dir, 'Pt4_tumour_growth_status.rds'))
  tg_535 <- readRDS(paste0(output_dir, 'Pt5_tumour_growth_status.rds'))
  tg_1035 <- readRDS(paste0(output_dir, 'Pt6_tumour_growth_status.rds'))
  
  pUnRx <- cowplot::plot_grid(collapse_trees$SA501,collapse_trees$SA530, collapse_trees$SA604, ncol=3, rel_widths = c(1,0.8,1),
                              labels = c('Pt1','Pt2','Pt3'), label_y = 0.2,label_x=0.8, hjust = -0.5, label_size=12)
  pSeries <- cowplot::plot_grid(collapse_trees$SA609,collapse_trees$SA535, collapse_trees$SA1035, ncol=3, rel_widths = c(0.8,1,1),
                                labels = c('Pt4','Pt5','Pt6'), label_y = 0.2,label_x=0.75, hjust = -0.5, label_size=12)
  
  # pRxDLP <- cowplot::plot_grid(NULL, res_barDLP_609$p, res_barDLP_535$p, res_barDLP_1035$p, ncol=4, rel_widths = c(0.05,1.15,1,1),
  #                              hjust = 0, label_size=12, labels = c('c',NULL, NULL, NULL)) + 
  #   theme(plot.background = element_rect(fill = "white", colour = "white"))
  tg_609 <- cowplot::plot_grid(NULL, tg_609, NULL, nrow=1, rel_widths = c(0.04,1,0.04))
  tg_535 <- cowplot::plot_grid(NULL, tg_535, NULL, nrow=1, rel_widths = c(0.04,1,0.04))
  tg_1035 <- cowplot::plot_grid(NULL, tg_1035, NULL, nrow=1, rel_widths = c(0.04,1,0.04))
  p4 <- cowplot::plot_grid(res_barDLP_609$p, tg_609, nrow=2, rel_heights = c(2,0.8))
  p5 <- cowplot::plot_grid(res_barDLP_535$p, tg_535, nrow=2, rel_heights = c(2,0.8))
  p6 <- cowplot::plot_grid(res_barDLP_1035$p, tg_1035, nrow=2, rel_heights = c(2,0.8))
  
  pRxDLP <- cowplot::plot_grid(NULL,p4, p5, p6, ncol=4, rel_widths = c(0.05,1.15,1,1),
                               hjust = 0, label_size=12, labels = c('c',NULL, NULL, NULL)) + 
    theme(plot.background = element_rect(fill = "white", colour = "white"))
                               
  # pRxDLP <- cowplot::plot_grid(NULL, res_barDLP_609$p, res_barDLP_535$p, res_barDLP_1035$p, ncol=4, rel_widths = c(0.05,1.15,1,1),
  #                              hjust = 0, label_size=12, labels = c('c',NULL, NULL, NULL)) + 
  #   theme(plot.background = element_rect(fill = "white", colour = "white"))

  
    
  # ptotal_DLP <- cowplot::plot_grid(NULL, pUnRx,NULL, pSeries, pRxDLP, nrow=5, rel_heights = c(0.1,1,0.1,1,2),
  #                                  labels = c('     -Rx PDX tumors','',
  #                                             '     -Rx/+Rx time series PDX tumors',''),
  #                                  # , label_x = 0, label_y = 1
  #                                  label_x = .01, hjust = 0, label_size=13)+
  #   theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  pDLP_tree1 <- cowplot::plot_grid(NULL, pUnRx,NULL, pSeries, nrow=4, rel_heights = c(0.1,1,0.1,1),
                                   labels = c('     Hierarchical dynamics of clones in Untreated PDX','',
                                              '     Hierarchical dynamics of clones in PDX time series',''),
                                   #, label_x = 0, label_y = 1
                                   label_x = .01, hjust = 0, label_size=11)+
    # c('     UnRx PDX tumors','',
    #   '     UnRx/Rx,RxH time series PDX tumors','')
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  output_dir <- "/home/htran/Projects/farhia_project/drug_resistant_material/materials/umap_figs/main_fig2/"
  p_manhattan <- readRDS(paste0(output_dir,"median_CN_distance_6series.rds"))
  
  
  
  pDLP_tree <- cowplot::plot_grid(pDLP_tree1, p_manhattan, nrow=1, rel_widths = c(1,0.2), hjust = 0, label_size=12, labels = c('a',' b')) + #p_manhattan
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  # ptotal_DLP <- cowplot::plot_grid(pDLP_tree, pRxDLP, nrow=2,rel_heights = c(2.2,2))
  # ptotal_DLP
  
  
  
}

## All plots here
plot_clonealign_correlation <- function(){
  # save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/'
  # base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/'
  output_dir <- "/home/htran/Projects/farhia_project/drug_resistant_material/materials/umap_figs/main_fig2/"
  p_corr <- readRDS(paste0(output_dir,"clonealign_correlation-SA609.rds"))
  res <- readRDS(paste0(output_dir,"clonealign_correlation-SA609_legends.rds"))
  # p1 <- p1 + labs(title='Pt4 - clonal proportions') + theme(axis.title = element_text(size=9, face="bold",family=my_font, hjust = 0.5))
  # pSA609_clonealign <- cowplot::plot_grid(p1, res$colplt,res$sizeplt, ncol = 1, rel_heights = c(10,0.95,0.9))
  
  # base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/'
  # save_dir <- base_dir
  # 
  # Clonealign Correlation for all patients
  # save_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/"
  # stat <- data.table::fread(paste0(save_dir, '10x_dlp_stat_correlation.csv')) %>% as.data.frame()
  # stat$Pt <- stat$PDX
  # my_font <- "Helvetica"
  # p_corr_all <- ggplot(stat, aes(x = DLP, y = clonealign)) + 
  #   geom_point(aes(color = Pt, shape=Pt), alpha=0.9, size=3)  #, size=2*log10(pct_genes)size=4.5
  # 
  # p_corr_all <- p_corr_all + thesis_theme
  # lg_pos <- "top"
  # # lg_pos <- c(0.7, 0.2)
  # p_corr_all <- p_corr_all + ggplot2::theme(legend.position = lg_pos,
  #                         legend.title = element_blank(), 
  #                         legend.key=element_blank(),
  #                         legend.text = element_text(color="black",size=9, hjust = 0.5, family=my_font)) +
  #   labs(x='Clonal proportion DLP', y='Clonal proportion clonealign 10x', 
  #        title = '6 patients - Clonal proportions')
  # 
  # p_corr_all <- p_corr_all + guides(color = guide_legend(nrow = 2, override.aes = list(size=4)))
  # saveRDS(p_corr_all, paste0(save_dir,'clonealign_correlation_all_patients.rds'))
  # save_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/clonealign_plot/"
  p_corr_all <- readRDS(paste0(output_dir,'clonealign_correlation_all_patients.rds'))
  # stat_DLP <- cowplot::plot_grid(res$colplt,res$colplt, ncol = 1, rel_heights = c(2, 1))
  pSA609_clonealign <- cowplot::plot_grid(p_corr_all,  res$colplt, res$sizeplt + 
                                            theme(legend.title = element_blank())
                                          , ncol = 1, rel_heights = c(4,0.4,0.4))
  fig2_part2 <- cowplot::plot_grid(pSA609_clonealign, p_corr, ncol = 2, rel_widths = c(1,1), 
                                   hjust = 0, label_size=12, labels = c('e','f'))+
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  
  ## TO DO add legend here 
  plg_tg <- readRDS(paste0(output_dir, 'tumour_growth_legend.rds'))
  p_fig2 <- cowplot::plot_grid(pDLP_tree, pRxDLP, plg_tg, fig2_part2, ncol = 1, 
                               rel_heights = c(0.6,0.6,0.1,0.8))+
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/materials/'
  output_dir <- paste0(base_dir,'umap_figs/main_fig2/')
  # dir.create(output_dir)
  ggsave(paste0(output_dir,"Fig2_DLPtree_clonealign_correlation.png"),
         plot = p_fig2,
         height = 11,
         width = 8,
         type = "cairo-png",
         dpi=150
  )
  # BiocManager::install('svglite')
  ggsave(paste0(output_dir,"Fig2_DLPtree_clonealign_correlation_v2.svg"),
         plot = p_fig2,
         height = 12,
         width = 8,
         # type = "cairo-png",
         dpi=250
  )
  # ggsave(paste0(output_dir,"Fig2_DLPtree_clonealign_correlation.pdf"),
  #        plot = p_fig2,
  #        height = 11,
  #        width = 7.5,
  #        useDingbats=F,
  #        dpi = 150
  # )
  # 
  
  
}







# old version
plotDLP_barplot_v1 <- function(){
  # c('UnRx','Rx','RxH')
  tag <- 'SA609'
  library_grouping_fn <- paste0(base_dir,'cell_clones/',tag, '_DLP_library_groupings.csv')
  cell_clones_fn <- paste0(base_dir,'cell_clones/',tag, '_cell_clones_total.csv')
  cell_clones <- get_meta_data(cell_clones_fn, library_grouping_fn, tag)
  # data.table::fwrite(cell_clones, paste0(base_dir,'cell_clones/',tag, '_cell_clones_metadata.csv'))
  
  # Look at clone E
  df <- color_df %>%
    dplyr::filter(datatag==tag)
  cols_use <- df$colour
  names(cols_use) <- df$clone_id
  unique(cell_clones$treatment_desc)
  res_SA609 <- plot_fill_barplot_wholedataset(cell_clones, cols_use, output_dir, 
                                              datatag, plottitle=NULL, plotlegend=F)
  res_SA609$p
  res_SA609$plg
  
  
  tag <- 'SA535'
  library_grouping_fn <- paste0(base_dir,'cell_clones/',tag, '_DLP_library_groupings.csv')
  cell_clones_fn <- paste0(base_dir,'cell_clones/',tag, '_cell_clones.csv')
  cell_clones <- get_meta_data(cell_clones_fn, library_grouping_fn, tag)
  # data.table::fwrite(cell_clones, paste0(base_dir,'cell_clones/',tag, '_cell_clones_metadata.csv'))
  df <- color_df %>%
    dplyr::filter(datatag==tag)
  cols_use <- df$colour
  names(cols_use) <- df$clone_id
  res_SA535 <- plot_fill_barplot_wholedataset(cell_clones, cols_use, output_dir, 
                                              datatag, plottitle=NULL, plotlegend=F)
  
  # res_SA535$p
  # res_SA535$plg
  
  tag <- 'SA1035'
  library_grouping_fn <- paste0(base_dir,'cell_clones/',tag, '_DLP_library_groupings.csv')
  cell_clones_fn <- paste0(base_dir,'cell_clones/',tag, '_cell_clones.csv')
  cell_clones <- get_meta_data(cell_clones_fn, library_grouping_fn, tag)
  # data.table::fwrite(cell_clones, paste0(base_dir,'cell_clones/',tag, '_cell_clones_metadata.csv'))
  df <- color_df %>%
    dplyr::filter(datatag==tag)
  cols_use <- df$colour
  names(cols_use) <- df$clone_id
  res_SA1035 <- plot_fill_barplot_wholedataset(cell_clones, cols_use, output_dir, 
                                               datatag, plottitle=NULL, plotlegend=F)
  # res_SA1035$p
  # res_SA1035$plg
  
}