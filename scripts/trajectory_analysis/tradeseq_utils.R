suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot, quietly = TRUE)
  # library(mclust, quietly = TRUE)
  library(tradeSeq)
  # library(grDevices)
  library(RColorBrewer)
  # library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(BiocParallel)
  library(tidyverse)
})
suppressPackageStartupMessages(library(circlize))
options(dplyr.summarise.inform = FALSE)
options(tidyverse.quiet = TRUE)
my_font <- "Helvetica"

get_pretty_gene_type_labels <- function(obs_genes_df, input_dir=NULL){
  if(is.null(input_dir)){
    input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/trajectory_genes/'  
  }
  meta_gm <- data.table::fread(paste0(input_dir,'meta_gene_module_labels_',datatag,'.csv.gz'))
  # if(datatag=='SA609'){
  #   
  #   
  # }else if(datatag=='SA535'){
  #   
  # }else{
  #   
  # }
  obs_genes_df <- obs_genes_df %>%
    dplyr::left_join(meta_gm, by=c('gene_type'='gene_type_module')) %>%
    dplyr::rename(gene_type_origin=gene_type, gene_type=gm_manuscript_lb)
  print(summary(as.factor(obs_genes_df$gene_type)))
  return(obs_genes_df)
  
}

## palette_use can be "Set2", "Accent", "Dark2", ...
## colorRampPalette if you want to use more than the number of colors in the palette
get_color_code <- function(nb.cols, palette_use="Dark2"){
  
  # brewer.pal(8, "Accent")
  # display.brewer.pal(8, "Accent")
  # display.brewer.pal(8, "Dark2")
  # mycolors <- colorRampPalette(brewer.pal(8, palette_use))(nb.cols)
  mycolors <- brewer.pal(8, palette_use)
  mycolors <- mycolors[1:nb.cols]
  return(mycolors)
}
get_lineage_colors_linetype <- function(df, datatag){
  if(datatag=='SA609'){
    df <- df %>%
      dplyr::mutate(lineage=case_when(
        lineage=='L1' ~ "L3-Rx",
        lineage=='L2' ~ "L2-RxH",
        lineage=='L3' ~ "L1-UnRx",
        TRUE ~ lineage
      ))
    
    lineage_ls <- c('L1-UnRx','L2-RxH','L3-Rx')
    # cols_use <- c('#717171','#ffd500','#0000FF') # darkgrey, #9acd32=yellowgreen, darkblue
    cols_use <- get_color_code(length(lineage_ls), palette_use="Dark2")
    # names(cols_use) <- lineage_ls
    linetypes <- c('solid','longdash','twodash')
    # names(linetypes) <- lineage_ls
    line_sizes <- c(0.8, 1.2, 1)
    # names(line_sizes) <- lineage_ls
  }else if(datatag=='SA535'){
    df <- df %>%
      dplyr::mutate(lineage=case_when(
        lineage=='L1' ~ "L3-Rx,RxH",
        lineage=='L2' ~ "L2-Rx,RxH",
        lineage=='L4' ~ "L1-UnRx",
        lineage=='L3' ~ "L4-1Rx",
        TRUE ~ lineage
      ))
    
    
    lineage_ls <- c('L1-UnRx','L2-Rx,RxH','L3-Rx,RxH','L4-1Rx')
    # cols_use <- c('#717171','#9ACD32','#402994','#00FFFF') # darkgrey, #9acd32=yellowgreen, darkblue, cyan
    cols_use <- get_color_code(length(lineage_ls), palette_use="Dark2")
    # names(cols_use) <- lineage_ls
    linetypes <- c('solid','longdash','twodash','dotted')
    # names(linetypes) <- lineage_ls
    line_sizes <- c(0.9, 1.2, 1, 0.7)
    # names(line_sizes) <- lineage_ls
  }else if(datatag=='SA1035'){  ## SA1035
    # print('')
    df <- df %>%
      dplyr::mutate(lineage=case_when(
        lineage=='L1' ~ "L4-Rx,RxH",
        lineage=='L2' ~ "L3-UnRx,RxH",
        lineage=='L4' ~ "L2-Rx,RxH",
        lineage=='L3' ~ "L1-1Rx",
        TRUE ~ lineage
      ))
    
    lineage_ls <- c('L3-UnRx,RxH','L4-Rx,RxH','L2-Rx,RxH','L1-1Rx')
    # cols_use <- c('#717171','yellowgreen','#402994','cyan') # darkgrey, #9acd32=yellowgreen, darkblue
    cols_use <- get_color_code(length(lineage_ls), palette_use="Dark2")
    # names(cols_use) <- lineage_ls
    linetypes <- c('solid','longdash','twodash','dotted')
    # names(linetypes) <- lineage_ls
    line_sizes <- c(0.9, 1.3, 1, 0.7)
    # names(line_sizes) <- lineage_ls
  }else{
    print('Pls set color, line type for this series')
    return(NULL)
  }
  plt_setting <- tibble::tibble(lineage=lineage_ls, used_color=cols_use,
                                linetypes=linetypes, line_sizes=line_sizes)
  return(list(df=df, plt_setting=plt_setting))
  
}
# , xlimit=c(0,50)
viz_given_smooth_exp <- function(df, cols_use, line_sizes, x_text_disp='none', x_tick_disp='none',
                                 xplt='pseudotime',yplt='gene_exp', colorplt='lineage'){
  my_font <- "Helvetica" 
  if(x_text_disp=='none'){
    disp <- element_blank()
  }else{
    disp <- element_text(size=13, family = my_font, color='black')
  }
  
  
  p <- ggplot2::ggplot(df, aes_string(x=xplt, y=yplt)) +#, linetype=colorplt 
    # ggplot2::geom_point(alpha = 0.15, aes_string(colour = colorplt)) + 
    # ggplot2::geom_density_2d(alpha = 0.4, aes_string(colour = colorplt)) +#
    # ggplot2::geom_density_2d_filled(contour_var = "density", 
    #                                 alpha = 0.5, aes_string(colour = colorplt)) +
    # stat_density_2d(
    #   geom = "raster",
    #   aes_string(colour = colorplt),
    #   contour = FALSE
    # ) +
    ggplot2::geom_smooth(span = 0.8, method = 'loess', se = F, 
                         aes_string(colour = colorplt, size=colorplt)) + #, size=colorplt, data=stat_df, aes_string(x=xplt, y=yplt, colour=colorplt). method='lm', size=1.3, formula = y ~ poly(x, 3)
    scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) + 
    theme_bw() + 
    scale_color_manual(values = cols_use)+
    # scale_linetype_manual(values=linetypes)+
    scale_size_manual(values=line_sizes)+
    # xlim(0, 50) + 
    theme(#strip.text.y = element_text(angle = 0),
      # legend.position = 'bottom',
      legend.position = 'none',
      legend.text=element_text(size=11, family = my_font, color='black'),
      legend.title=element_text(size=11, family = my_font, color='black'),
      # axis.title = element_text(size=12, family = my_font, color='black'),
      axis.title = element_blank(),
      # axis.text.y = element_text(size=6, family = my_font, color='black'),
      axis.text.y = element_blank(),
      axis.text.x = disp,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "null"),
      panel.spacing = unit(c(0, 0, 0, 0), "null"),
      panel.background = element_blank(),
      strip.background = element_rect(color='white', fill='white'),
      # strip.text.y = element_text(size=12, family = my_font, color='black'),
      strip.text.y = element_blank(),
      axis.line.x.bottom = element_line(color='darkgray', size = 0.3),
      axis.line.y.left = element_line(color='darkgray', size = 0.3)) + 
    labs(y='Relative expression', x = 'Pseudotime', title=NULL)
  if(x_tick_disp=='none'){
    p <- p + theme(axis.ticks.x = element_blank())
  }
  return(p)
}
# output_dir: contains files: mtx_hm.csv.gz, obs_genes_hm.csv.gz, obs_cells_hm.csv.gz
viz_gene_modules_density_plot <- function(xplt='pseudotime',yplt='gene_exp',
                                          colorplt='lineage',datatag='',
                                          input_dir, output_dir=NULL){
 
  if(is.null(output_dir)){
    output_dir <- paste0(input_dir,'figs/')
  }
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  # save_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/"
  # ## phm <- viz_heatmap(ts_sce, genes_df, save_dir, datatag, plttitle)
  # output_dir <- save_dir
  
  
  ## SA535
  # output_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/'
  # datatag <- 'SA535'
  
  # save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/')
  # ts_sce <- readRDS(paste0(save_dir,'tradeSeq_3000/', "fitGAM_out.rds"))
  # nb_estimated_points <- 100
  # yhatSmooth <- tradeSeq::predictSmooth(ts_sce, #[,rownames(metacells)]
  #                                       gene = obs_genes_df$ens_gene_id, 
  #                                       nPoints = nb_estimated_points)
  # # class(yhatSmooth)
  # # unique(yhatSmooth$lineage)
  # obs_lineages <- NULL
  # if(!is.null(obs_lineages)){  ##ex: obs_lineages <- c(1,2,3,4)
  #   yhatSmooth <- yhatSmooth %>%
  #     dplyr::filter(lineage %in% obs_lineages)
  # }
  # yhatSmooth$desc <- paste0(yhatSmooth$lineage,'_', yhatSmooth$time)
  
  
  exp_mtx <- data.table::fread(paste0(input_dir,"mtx_hm.csv.gz")) %>% as.data.frame()
  obs_genes_df <- data.table::fread(paste0(input_dir,"obs_genes_hm.csv.gz")) %>% as.data.frame()
  obs_cells_df <- data.table::fread(paste0(input_dir,"obs_cells_hm.csv.gz")) %>% as.data.frame()
  # head(obs_cells_df)
  rownames(exp_mtx) <- exp_mtx$ens_gene_id
  # unique(obs_genes_df$gene_type)
  # exp_mtx$ens_gene_id <- NULL
  # dim(exp_mtx)
  # df <- exp_mtx[gm$ens_gene_id,1:100]
  # df <- exp_mtx
  # unique(obs_genes_df$gene_type)
  obs_lineages <- rep(1:length(unique(obs_cells_df$lineage))) 
  # if(datatag=='SA609'){
  #   obs_lineages <- rep(1:length(unique(obs_cells_df$lineage))) # SA609  
  # }else if(datatag=='SA535'){
  #   # obs_lineages <- c(1,2,4) # SA535
  # }else{
  #   obs_lineages <- rep(1:length(unique(obs_cells_df$lineage)))
  # }
  
  
  # gtm_ls <- gsub('Module','',unique(obs_genes_df$gene_type))
  
  df <- tibble::tibble()
  for(l in obs_lineages){
    cols_use <- colnames(exp_mtx)[grepl(paste0(as.character(l),'_'),colnames(exp_mtx))]
    cols_use <- c('ens_gene_id',cols_use)
    tmp <- exp_mtx %>%
      dplyr::select(all_of(cols_use)) %>%
      tidyr::pivot_longer(!ens_gene_id, names_to = 'pseudotime', values_to = 'gene_exp')#cols = colnames(df),
    tmp$lineage <- paste0('L',l)
    df <- dplyr::bind_rows(df, tmp)
  }
  # dim(df)
  # colnames(obs_genes_df)
  # obs_genes_df$description[1]
  obs_genes_df <- get_pretty_gene_type_labels(obs_genes_df)
  obs_genes_df$gene_type <- gsub('Module','M',obs_genes_df$gene_type)
  genes_df <- obs_genes_df %>%
    dplyr::select(ens_gene_id,gene_type)
  df <- df %>% left_join(genes_df, by='ens_gene_id')  
  # unique(df$lineage)
  # if(datatag=='SA1035'){
  #   df <- df %>%
  #     dplyr::filter(lineage!='L5')
  # }
  res <- get_lineage_colors_linetype(df, datatag)
  df <- res$df
  cols_use <- res$plt_setting$used_color
  names(cols_use) <- res$plt_setting$lineage
  color_df <- tibble::tibble(color_lineage=cols_use, lineage=names(cols_use))
  data.table::fwrite(color_df, paste0(input_dir, datatag,'_smooth_line_lineage_colors.csv.gz'))
  
  
  linetypes <- res$plt_setting$linetypes
  names(linetypes) <- res$plt_setting$lineage
  line_sizes <- res$plt_setting$line_sizes
  names(line_sizes) <- res$plt_setting$lineage
  ##SA535
  # meta_lineage <- data.frame(lineage=c("Lineage 1","Lineage 2","Lineage 3","Lineage 4"),
  #                            lineage_desc=c("L3-Rx,RxH","L2-Rx,RxH",
  #                                           "L4-Rx,1Rx","L1-UnRx"))
  
  # df <- df %>%
  #   tidyr::pivot_longer(!ens_gene_id, names_to = 'pseudotime', values_to = 'gene_exp')#cols = colnames(df), 
  # head(df)
  # dim(df)
  
  
  # df$pseudotime <- 
  ps <- lapply(strsplit(df$pseudotime,'_'), function(x){
    return(x[2])    
  })
  df$pseudotime <- as.numeric(unlist(ps))
  # df$pseudotime[1]
  
  # stat_df <- tibble::tibble()
  # gms <- unique(df[,colorplt])
  # for(t in gms){
  #   tmp <- df %>%
  #     dplyr::filter(!!sym(colorplt)==t)%>%
  #     # dplyr::filter(x < as.numeric(quantile(x, probs = 0.95)) & x>as.numeric(quantile(x, probs = 0.05)))
  #     dplyr::mutate(gene_exp=case_when(
  #       (gene_exp > as.numeric(quantile(gene_exp, probs = 0.95))) ~ as.numeric(quantile(gene_exp, probs = 0.95)),
  #       (gene_exp < as.numeric(quantile(gene_exp, probs = 0.05))) ~ as.numeric(quantile(gene_exp, probs = 0.05)),
  #       TRUE ~ gene_exp
  #     ))
  #   # tmp <- df %>%
  #   #   dplyr::filter(!!sym(colorplt)==t)%>%
  #   #   # dplyr::filter(x < as.numeric(quantile(x, probs = 0.95)) & x>as.numeric(quantile(x, probs = 0.05)))
  #   #   dplyr::mutate(!!sym(xplt)=case_when(
  #   #     (!!sym(xplt) > as.numeric(quantile(!!sym(xplt), probs = 0.95))) ~ as.numeric(quantile(!!sym(xplt), probs = 0.95)),
  #   #     (!!sym(xplt) < as.numeric(quantile(!!sym(xplt), probs = 0.05))) ~ as.numeric(quantile(!!sym(xplt), probs = 0.05)),
  #   #     TRUE ~ !!sym(xplt)
  #   #   ))
  #   stat_df <- dplyr::bind_rows(stat_df, tmp)
  # }
  # dim(stat_df)
  
  # stat_df <- stat_df %>% 
  #   dplyr::mutate(bin = cut_interval(!!sym(xplt), n = 200))%>% 
  #   dplyr::group_by(!!sym(colorplt), bin)%>%  
  #   dplyr::summarise(!!sym(yplt)=median(!!sym(yplt)), !!sym(xplt) = median(!!sym(xplt)))%>% 
  #   dplyr::ungroup()%>% 
  #   dplyr::select(-bin)
  # stat_df <- df
  # exp_095 <- as.numeric(quantile(gene_exp, probs = 0.95))
  # exp_005 <- as.numeric(quantile(gene_exp, probs = 0.05))
  # stat_df <- stat_df %>%
  #   #     # dplyr::filter(x < as.numeric(quantile(x, probs = 0.95)) & x>as.numeric(quantile(x, probs = 0.05)))
  #       dplyr::mutate(gene_exp=case_when(
  #         (gene_exp > exp_095) ~ exp_095,
  #         (gene_exp < exp_005) ~ exp_005,
  #         TRUE ~ gene_exp
  #       ))
  # stat_df <- stat_df %>% 
  #   # dplyr::mutate(bin = cut_interval(pseudotime, n = 200))%>% 
  #   dplyr::group_by(lineage, pseudotime, gene_type)%>%  
  #   dplyr::summarise(pseudotime=mean(pseudotime), gene_exp = mean(gene_exp))%>% 
  #   dplyr::ungroup()#%>% 
  #   # dplyr::select(-bin)
  
  # Show the contour only
  # unique(df$lineage)
  dim(df)
  summary(as.factor(df$gene_type))
  xmax <- round(max(df$pseudotime)/10) *10
  plt_ls <- list()
  gene_modules <- unique(df$gene_type)
  for(gt in gene_modules){
    tmp <- df %>% 
      dplyr::filter(gene_type==gt)
    plt_ls[[gt]] <- viz_given_smooth_exp(tmp, cols_use, line_sizes, x_text_disp='none', x_tick_disp='none', #xlimit=c(0,xmax),
                                         xplt, yplt, colorplt)
  }
  # drawing x axis values
  gt <- paste0('M',length(gene_modules))
  tmp <- df %>% 
    dplyr::filter(gene_type==gt)
  plt_ls[[gt]] <- viz_given_smooth_exp(tmp, cols_use, line_sizes, x_text_disp='display', x_tick_disp='display', #xlimit=c(0,xmax),
                                       xplt, yplt, colorplt)
  # names(plt_ls)
  plt_ls2 <- list()
  for(m in rep(1:length(gene_modules))){
    plt_ls2[[paste0('M',m)]] <- plt_ls[[paste0('M',m)]]
  }
  
  p <- cowplot::plot_grid(plotlist = plt_ls2, ncol=1, align='v') #, labels = names(plt_ls2)
  
  # p <- ggplot2::ggplot(df, aes_string(x=xplt, y=yplt)) +#, linetype=colorplt , size=colorplt
  #   # ggplot2::geom_point(alpha = 0.15, aes_string(colour = colorplt)) + 
  #   # ggplot2::geom_density_2d(alpha = 0.4, aes_string(colour = colorplt)) +#
  #   # ggplot2::geom_density_2d_filled(contour_var = "density", 
  #   #                                 alpha = 0.5, aes_string(colour = colorplt)) +
  #   # stat_density_2d(
  #   #   geom = "raster",
  #   #   aes_string(colour = colorplt),
  #   #   contour = FALSE
  #   # ) +
  #   ggplot2::geom_smooth(span = 0.8, method = 'loess', se = F, 
  #                        aes_string(colour = colorplt)) + #, size=colorplt, data=stat_df, aes_string(x=xplt, y=yplt, colour=colorplt). method='lm', size=1.3, formula = y ~ poly(x, 3)
  #   scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) + 
  #   # facet_grid(rows = vars(gene_type)) + #, scales = "free", space = "free"
  #   facet_wrap(vars(gene_type), nrow=length(unique(df$gene_type))) +  #, scales = "free_x"
  #     theme_bw() + 
  #     scale_color_manual(values = cols_use)+
  #     # scale_linetype_manual(values=linetypes)+
  #     # scale_size_manual(values=line_sizes)+
  #     theme(#strip.text.y = element_text(angle = 0),
  #       legend.position = 'bottom',
  #       legend.text=element_text(size=11, family = my_font, color='black'),
  #       legend.title=element_text(size=11, family = my_font, color='black'),
  #       axis.title = element_text(size=12, family = my_font, color='black'),
  #       # axis.text.y = element_text(size=6, family = my_font, color='black'),
  #       axis.text.y = element_blank(),
  #       axis.text.x = element_text(size=11, family = my_font, color='black'),
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(),
  #       panel.border = element_blank(),
  #       plot.margin = unit(c(0, 0, 0, 0), "null"),
  #       panel.spacing = unit(c(0, 0, 0, 0), "null"),
  #       panel.background = element_blank(),
  #       strip.background = element_rect(color='white', fill='white'),
  #       # strip.text.y = element_text(size=12, family = my_font, color='black'),
  #       strip.text.y = element_blank(),
  #       axis.line.x.bottom = element_line(color='darkgray', size = 0.3),
  #       axis.line.y.left = element_line(color='darkgray', size = 0.3)) + 
  #   labs(y='Relative expression', x = 'Pseudotime', title=NULL)#x = 'Pseudotime'
  # # p
  # # p
  # lg <- cowplot::get_legend(p)
  # plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  # p <- p + theme(legend.position = 'none')
  # png(paste0(output_dir,"gene_module_exp_pseudotime_withlegend.png"),
  #     height = 2*700, width=2*310, res = 2*72)
  # print(p)
  # dev.off()
  # png(paste0(output_dir,"gene_module_exp_pseudotime_withoutlegend.png"),
  #     height = 2*700, width=2*250, res = 2*72)
  # print(p)
  # dev.off()
  # png(paste0(output_dir,"density_plots.png"),height = 2*300, width=2*800,res = 2*72)
  # print(p)
  # dev.off()
  # ggsave(paste0(output_dir,"density_plot_",datatag, ".png"),
  #        plot = p,
  #        height = 5,
  #        width = 2.2,
  #        # useDingbats=F,
  #        dpi=150)
  # ggsave(paste0(output_dir,"density_plot_",datatag, "_lg.png"),
  #        plot = plg,
  #        height = 0.6,
  #        width = 4,
  #        # useDingbats=F,
  #        dpi=150)
  saveRDS(p, paste0(output_dir,"density_plot_",datatag, ".rds"))
  plg <- NULL
  # saveRDS(plg, paste0(output_dir,"density_plot_",datatag, "_lg.rds"))
  
  
  return(list(p_smooth_lineage_exp=p, plg=plg))
}
# ts_sce: tradeSeq obj
# obs_genes_df: ens_gene_id, gene_symbol, gene_type (ex: Module1, Module2,..)
# tag <- 'clusters_activated_repressed'
viz_heatmap <- function(ts_sce, obs_genes_df, output_dir, datatag, plttitle, obs_lineages=NULL){
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  meta_genes <- get_meta_genes(obs_genes_df)
  nb_estimated_points <- 100
  yhatSmooth <- tradeSeq::predictSmooth(ts_sce, #[,rownames(metacells)]
                                        gene = obs_genes_df$ens_gene_id, 
                                        nPoints = nb_estimated_points)
  # class(yhatSmooth)
  # unique(yhatSmooth$lineage)
  if(!is.null(obs_lineages)){  ##ex: obs_lineages <- c(1,2,3,4)
    yhatSmooth <- yhatSmooth %>%
      dplyr::filter(lineage %in% obs_lineages)
  }
  yhatSmooth$desc <- paste0(yhatSmooth$lineage,'_', yhatSmooth$time)
  l1 <- yhatSmooth %>%
    # dplyr::arrange(lineage)%>%
    dplyr::select(-lineage, -time)%>%
    tidyr::pivot_wider(names_from = 'desc', values_from = 'yhat')%>%
    tibble::column_to_rownames('gene')
  # dim(l1)
  lgs <- unique(yhatSmooth$lineage)
  lineages <- c()
  for(lg in lgs){
    lineages <- c(lineages, rep(paste0('Lineage ',lg), nb_estimated_points))
  }
  obs_cells_df <- data.frame(lineage=lineages)
  
  meta_lineage <- get_meta_lineages_desc(datatag)
  if(!is.null(meta_lineage)){
    obs_cells_df <- obs_cells_df %>%
      inner_join(meta_lineage, by='lineage')%>%
      select(-lineage)%>%
      dplyr::rename(lineage=lineage_desc)  
  }
  
  ## Chop the max, and min gene expression to this max value --> heatmap colors
  max_exp <- 4
  exp_mtx <- t(scale(t(l1)))
  exp_mtx[exp_mtx>max_exp] <- max_exp
  exp_mtx[exp_mtx<-max_exp] <- -max_exp
  # rownames(exp_mtx) <- meta_genes[rownames(exp_mtx),'gene_symbol']
  dim(exp_mtx)
  dim(obs_cells_df)
  rownames(obs_cells_df) <- colnames(l1)
  
  
  p <- viz_genes_exp_lineages_hm(exp_mtx, obs_genes_df, obs_cells_df, output_dir, plttitle)
  exp_mtx <- as.data.frame(exp_mtx)
  exp_mtx$ens_gene_id <- rownames(exp_mtx)
  # sum(rownames(exp_mtx)==obs_genes_df$ens_gene_id)
  data.table::fwrite(exp_mtx, paste0(output_dir,"mtx_hm.csv.gz"))
  data.table::fwrite(obs_genes_df, paste0(output_dir,"obs_genes_hm.csv.gz"))
  obs_cells_df$cell_id <- rownames(obs_cells_df)
  data.table::fwrite(obs_cells_df, paste0(output_dir,"obs_cells_hm.csv.gz")) 
  
  return(p)
}

plot_legends_heatmap <- function(output_dir){
  library(ComplexHeatmap)
  library(circlize)
  # gt <- c('cis','trans','unmapped')
  # col_gt <- c("cis" = "chocolate", "trans" = "blue2", "unMapped"="lightgrey")
  gt <- c('cis','trans')
  col_gt <- c("cis" = "chocolate", "trans" = "blue2")
  lg_cistrans = ComplexHeatmap::Legend(labels = gt, title = "Gene Type",labels_gp = gpar(fontsize = 13),
                                       legend_gp = grid::gpar(fill = col_gt))
  hm_plg_cistrans <- grid::grid.grabExpr(ComplexHeatmap::draw(lg_cistrans))
  saveRDS(hm_plg_cistrans, paste0(output_dir,"gene_type_cistrans_hm_lg.rds"))
  
  col_fun = colorRamp2(c(-4,0,4), c("blue", "white", "red"))
  lgd = ComplexHeatmap::Legend(col_fun = col_fun, title = "Avg Exp", at = c(-4,-2,0, 2, 4), direction = "horizontal")
  hm_plg_exp <- grid.grabExpr(ComplexHeatmap::draw(lgd))
  

  gt_chromatin <- c('Active','Bivalent','Repressed')
  col_chromatin <- c("Active" = "#DE3163", "Bivalent" = "#00FFFF", "Repressed"="#6495ED", "NotSignf"="lightgrey")
  col_chromatin <- col_chromatin[gt_chromatin]
  lg_chrm = ComplexHeatmap::Legend(labels = gt_chromatin, title = "Chromatin Status",labels_gp = gpar(fontsize = 13),
                                       legend_gp = grid::gpar(fill = col_chromatin))
  hm_plg_chromatin <- grid::grid.grabExpr(ComplexHeatmap::draw(lg_chrm))
  #ex: cowplot::plot_grid(hm_plg, ncol=1)
  saveRDS(hm_plg_chromatin, paste0(output_dir,"gene_type_hm_legend.rds"))
  
  res <- list(hm_plg_cistrans=hm_plg_cistrans, hm_plg_exp=hm_plg_exp, hm_plg_chromatin=hm_plg_chromatin)
  return(res)
  
}
viz_genes_exp_lineages_cistrans_anno_hm <- function(exp_mtx, obs_genes_df, obs_cells_df, 
                                                    cistrans_anno, meta_clone_lg,
                                                    output_dir,plttitle){
  rownames(obs_genes_df) <- NULL
  obs_genes_df <- obs_genes_df %>% 
    # as.data.frame()%>% 
    # select(-nb_genes)%>%
    select(gene_type, ens_gene_id, chromatin_st) #%>% 
    # tibble::column_to_rownames('ens_gene_id')
  # head(obs_genes_df)
  
  if('cell_id' %in% colnames(obs_cells_df)){ #is.null(rownames(obs_cells_df)) & 
    obs_cells_df <- obs_cells_df %>%
      tibble::column_to_rownames('cell_id')
  }
  obs_cells_df <- obs_cells_df %>%
    left_join(meta_clone_lg, by='lineage') %>%
    select(-lineage) %>%
    rename(lineage=lineage_desc)
  dim(obs_cells_df)
  head(obs_cells_df)
  ## Adding number of genes for each module
  meta_genes <- obs_genes_df %>% 
    dplyr::group_by(gene_type) %>% 
    dplyr::summarise(nb_genes=n())
  obs_genes_df <- obs_genes_df %>% inner_join(meta_genes, by='gene_type')
  obs_genes_df$gene_type <- paste0(obs_genes_df$gene_type,'\n(',obs_genes_df$nb_genes,')')
  obs_genes_df$nb_genes <- NULL
  obs_genes_df$gene_type <- gsub('Module','M',obs_genes_df$gene_type)
  obs_genes_df <- obs_genes_df[gtools::mixedorder(obs_genes_df$gene_type),, drop=F]
  head(obs_genes_df)
  # unique(obs_genes_df$gene_type)
  
  # unique(obs_genes_df$gene_type)
  # left_anno = ComplexHeatmap::rowAnnotation(Cls = factor(obs_genes_df$gene_type),
  #                                           annotation_legend_param = list(
  #                                             Cls = list(direction = "vertical")))
  left_anno <- NULL
  # srn <- ifelse(nb_genes<=50,T,F)
  srn <- F
  anno_col <- NULL
  
  ## draw legend 
  # col_fun = colorRamp2(c(-4,0,4), c("blue", "white", "red"))
  # lgd = ComplexHeatmap::Legend(col_fun = col_fun, title = "AvgExp", at = c(-4,-2,0, 2, 4), direction = "horizontal")
  # hm_plg <- grid.grabExpr(ComplexHeatmap::draw(lgd))
  library(ComplexHeatmap)
  gt <- c('cis','trans','unmapped')
  col_gt <- c("cis" = "chocolate", "trans" = "blue2", "unmapped"="lightgrey")
  lg_cistrans = ComplexHeatmap::Legend(labels = gt, title = "Gene Type",labels_gp = gpar(fontsize = 13),
  legend_gp = grid::gpar(fill = col_gt))
  hm_plg <- grid::grid.grabExpr(ComplexHeatmap::draw(lg_cistrans))
  #ex: cowplot::plot_grid(hm_plg, ncol=1)
  saveRDS(hm_plg, paste0(output_dir,"gene_type_hm_legend.rds"))
  
  meta_clone_lg$lineage_desc <- gsub('Lineage ','L',meta_clone_lg$lineage_desc)
  cistrans_anno_df <- cistrans_anno %>%
            inner_join(meta_clone_lg, by='lineage') %>%
            dplyr::select(-lineage) %>% #-clone
            tidyr::pivot_wider(names_from='lineage_desc', values_from = 'gene_type') %>%
            tibble::column_to_rownames('ens_gene_id')
  # cistrans_anno_df <- cistrans_anno_df[rownames(exp_mtx),]
  cistrans_anno_df <- cistrans_anno_df[obs_genes_df$ens_gene_id,]
  # anno_row <- ComplexHeatmap::rowAnnotation(
  #   GeneType = cbind(Lineage1 = sample(gt, dim(obs_genes_df)[1], replace = TRUE), 
  #                    Lineage2 = sample(gt, dim(obs_genes_df)[1], replace = TRUE),
  #                    Lineage3 = sample(gt, dim(obs_genes_df)[1], replace = TRUE)), 
  #   
  #   col = list(GeneType = col_gt), simple_anno_size = unit(0.6, "cm"), show_legend = F
  # )
  lg_ls <- gtools::mixedsort(colnames(cistrans_anno_df))
  cistrans_anno_df <- cistrans_anno_df %>%
    dplyr::select(all_of(lg_ls))
  
  gt_chromatin <- c('Active','Bivalent','Repressed')
  col_chromatin <- c("Active" = "#DE3163", "Bivalent" = "#00FFFF", "Repressed"="#6495ED", "NotSignf"="lightgrey")
  col_chromatin <- col_chromatin[unique(obs_genes_df$chromatin_st)]
  # chromatin_st <- c(rep('Active',t), rep('Bivalent',dim(obs_genes_df)[1]-t))
  chromatin_st <- obs_genes_df$chromatin_st
  anno_row <- ComplexHeatmap::rowAnnotation(
    GeneType = as.matrix(cistrans_anno_df), 
    Chromatin = chromatin_st, gap = unit(1, "points"),
    col = list(GeneType = col_gt, Chromatin=col_chromatin), 
    simple_anno_size = unit(0.45, "cm"), show_legend = F,
    annotation_name_rot = 90,
    annotation_name_gp= gpar(fontsize = 8)
  )
  
  # class(exp_mtx)[1:3, 1:3]
  exp_mtx <- as.matrix(exp_mtx[obs_genes_df$ens_gene_id,])
  rownames(obs_genes_df) <- NULL
  
  genes_clusters <- obs_genes_df$gene_type
  names(genes_clusters) <- obs_genes_df$ens_gene_id
  # cls <- ComplexHeatmap::cluster_within_group(t(exp_mtx), genes_clusters)
  # cls[1:3]
  p <- ComplexHeatmap::Heatmap(exp_mtx, na_col = "white",
                               show_column_names=F,
                               show_row_names = F,
                               # cluster_rows=cluster_within_group(t(exp_mtx), genes_clusters),
                               cluster_rows=T, 
                               cluster_columns=F,
                               name = "Avg Exp", 
                               # row_order = sort(rownames(test)),
                               row_split= genes_clusters,
                               # row_title_rot = 0,
                               row_gap = unit(1.2, "mm"),
                               column_split = obs_cells_df, 
                               # column_title = "ddd",
                               column_gap = unit(1.2, "mm"),
                               column_names_gp = grid::gpar(fontsize = 11),
                               column_title_gp = grid::gpar(fontsize = 11),
                               row_names_gp = grid::gpar(fontsize = 11, hjust=0.5),
                               row_title_gp = grid::gpar(fontsize = 11, hjust=0.5),
                               row_title_rot = 0,
                               row_names_rot = 0,
                               top_annotation=anno_col,
                               left_annotation = left_anno,
                               right_annotation = anno_row,
                               # bottom_annotation = yaxis_label,
                               # cell_fun = cell_func,
                               # row_dend_reorder=F,
                               show_column_dend = F,
                               show_row_dend = F,
                               show_heatmap_legend = F,
                               # row_title = plttitle, 
                               heatmap_legend_param = list(direction = "horizontal")
                               # column_title = "Pseudotime"
                               # column_title_side = "bottom"
  )
  # p
  # add title
  # nb_genes <- dim(obs_genes_df)[1]
  # hz <- ifelse(nb_genes<=50,(nb_genes*10+10),700)
  hz <- 700
  # tag <- 'clusters_activated_repressed'
  # title <- 'clusters_drug_holiday_genes'
  plttitle <- gsub(' ','_',plttitle)
  
  
  # png(paste0(output_dir,plttitle,".png"), height = 2*hz, width=2*700, res = 2*72)
  # # print(p)
  # ComplexHeatmap::draw(p, annotation_legend_side = "bottom",
  #                      heatmap_legend_side = "bottom",
  #                      padding = unit(c(4, 4, 2, 2), "mm"))
  # dev.off()
  
  saveRDS(p,paste0(output_dir,plttitle,".rds"))
  return(p)
}

viz_genes_exp_lineages_hm <- function(exp_mtx, obs_genes_df, obs_cells_df, 
                                      output_dir,plttitle){
  rownames(obs_genes_df) <- NULL
  obs_genes_df <- obs_genes_df %>% 
    # as.data.frame()%>% 
    # select(-nb_genes)%>%
    select(gene_type, ens_gene_id) %>% 
    tibble::column_to_rownames('ens_gene_id')
  # head(obs_genes_df)
  
  if('cell_id' %in% colnames(obs_cells_df)){ #is.null(rownames(obs_cells_df)) & 
    obs_cells_df <- obs_cells_df %>%
      tibble::column_to_rownames('cell_id')
  }
  
  ## Adding number of genes for each module
  meta_genes <- obs_genes_df %>% 
    dplyr::group_by(gene_type) %>% 
    dplyr::summarise(nb_genes=n())
  obs_genes_df <- obs_genes_df %>% inner_join(meta_genes, by='gene_type')
  obs_genes_df$gene_type <- paste0(obs_genes_df$gene_type,'(',obs_genes_df$nb_genes,')')
  obs_genes_df$nb_genes <- NULL
  # unique(obs_genes_df$gene_type)
  
  # unique(obs_genes_df$gene_type)
  # left_anno = ComplexHeatmap::rowAnnotation(Cls = factor(obs_genes_df$gene_type),
  #                                           annotation_legend_param = list(
  #                                             Cls = list(direction = "vertical")))
  left_anno <- NULL
  # srn <- ifelse(nb_genes<=50,T,F)
  srn <- F
  anno_col <- NULL
  # anno_col <- ComplexHeatmap::HeatmapAnnotation(treatment=factor(lineages))
  # Get genes modules
  # preprocess_mat <- as.matrix(counts(ts_sce[obs_genes_df$genes_use,]))
  # gene_module_df <- get_genes_modules(preprocess_mat=preprocess_mat, 
  #                                     resolution=0.2)
  # data.table::fwrite(gene_module_df, paste0(output_dir,'gene_modules_',datatag,'.csv'))
  p <- ComplexHeatmap::Heatmap(as.matrix(exp_mtx), na_col = "white",
                               show_column_names=F,
                               show_row_names = F,
                               cluster_rows=T,cluster_columns=F,
                               name = "Avg Exp", 
                               # row_order = sort(rownames(test)),
                               row_split= obs_genes_df,
                               # row_title_rot = 0,
                               row_gap = unit(2, "mm"),
                               column_split = obs_cells_df, 
                               # column_title = "ddd",
                               column_gap = unit(2, "mm"),
                               column_names_gp = grid::gpar(fontsize = 11),
                               column_title_gp = grid::gpar(fontsize = 14),
                               row_names_gp = grid::gpar(fontsize = 8),
                               row_title_gp = grid::gpar(fontsize = 9),
                               show_heatmap_legend = F,
                               top_annotation=anno_col,
                               left_annotation = left_anno,
                               # cell_fun = cell_func,
                               # row_dend_reorder=F,
                               show_column_dend = FALSE,
                               show_row_dend = FALSE,
                               # row_title = plttitle, 
                               heatmap_legend_param = list(direction = "horizontal")
                               # column_title = "Pseudotime", 
                               # column_title_side = "bottom"
  )
  
  # add title
  nb_genes <- dim(obs_genes_df)[1]
  hz <- ifelse(nb_genes<=50,(nb_genes*10+10),600)
  # tag <- 'clusters_activated_repressed'
  # title <- 'clusters_drug_holiday_genes'
  plttitle <- gsub(' ','_',plttitle)
  png(paste0(output_dir,plttitle,".png"), height = 2*hz, width=2*600, res = 2*72)
  # print(p)
  ComplexHeatmap::draw(p, annotation_legend_side = "bottom",
                       heatmap_legend_side = "bottom",
                       padding = unit(c(4, 4, 2, 2), "mm"))
  dev.off()
  # 
  # p
  saveRDS(p,paste0(output_dir,plttitle,".rds"))
  return(p)
}
get_meta_lineages_desc <- function(datatag){
  meta_lineage <- NULL
  if(datatag=='SA609'){
    # meta_lineage <- data.frame(lineage=c("Lineage 1","Lineage 2","Lineage 3"),
    #                            lineage_desc=c("Lineage 3: Rx","Lineage 2: RxH","Lineage 1: UnRx"))
    meta_lineage <- data.frame(lineage=c("Lineage 1","Lineage 2","Lineage 3"),
                               lineage_desc=c("L3-Rx","L2-RxH","L1-UnRx"))
  }else if(datatag=='SA535'){
    # meta_lineage <- data.frame(lineage=c("Lineage 1","Lineage 2","Lineage 3","Lineage 4"),
    #                            lineage_desc=c("Lineage 3: Rx,RxH","Lineage 2: Rx,RxH",
    #                                           "Lineage 4: Rx,1Rx","Lineage 1: UnRx"))
    
    # meta_lineage <- data.frame(lineage=c("Lineage 1","Lineage 2","Lineage 3","Lineage 4"),
    #                            lineage_desc=c("L3: Rx,RxH","L2: Rx,RxH",
    #                                           "L4: Rx,1Rx","L1: UnRx"))
    meta_lineage <- data.frame(lineage=c("Lineage 1","Lineage 2","Lineage 3","Lineage 4"),
                               lineage_desc=c("L3-Rx,RxH","L2-Rx,RxH",
                                              "L4-1Rx","L1-UnRx"))
    
  }else if(datatag=='SA1035'){
    # meta_lineage <- data.frame(lineage=c("Lineage 1","Lineage 2","Lineage 3","Lineage 4","Lineage 5"), 
    #                            lineage_desc=c("L4:Rx,RxH","L3:UnRx(Sen),RxH",
    #                                           "L2:Rx(Res),RxH","L1:1Rx,Rx", "L5:Rx,UnRx")) 
    meta_lineage <- data.frame(lineage=c("Lineage 1","Lineage 2","Lineage 3","Lineage 4","Lineage 5"), 
                               lineage_desc=c("L4-Rx,RxH","L3-UnRx(Sen),RxH",
                                              "L2-Rx(Res),RxH","L1-1Rx,Rx", "L5-Rx,UnRx")) 
  }else{
    print('Check meta lineage info!!!')
  }
  return(meta_lineage)
}

# save_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/tradeSeq_3000/SA609_tradeSeq/'
# output_fn <- paste0(save_dir,'sigf_gene_exp.csv.gz')
get_average_gene_exp_per_lineage <- function(datatag, ts_sce, obs_genes_ls, output_fn,
                                             nEstimatedPoints=100, save_data=T){
  meta_lineage <- get_meta_lineages_desc(datatag)
  sigf_gene_exp <-  tradeSeq::predictSmooth(models = ts_sce, 
                                            gene = obs_genes_ls, 
                                            nPoints = nEstimatedPoints)
  
  sigf_gene_exp$lineage <- paste0('Lineage ',sigf_gene_exp$lineage)
  # unique(sigf_gene_exp$lineage)
  
  sigf_gene_exp <- sigf_gene_exp %>% inner_join(meta_lineage, by=c('lineage'))
  # dim(sigf_gene_exp)
  # head(sigf_gene_exp)
  if(save_data){
    data.table::fwrite(sigf_gene_exp, output_fn)
    # sigf_gene_exp <- data.table::fread(output_fn) %>% as.data.frame()
  }
  return(sigf_gene_exp)
  
}

# avg_gene_exp <- get_average_gene_exp_per_lineage(datatag, ts_sce, obs_genes_ls, output_fn, nEstimatedPoints=100, save_data=T)
viz_given_gene_exp_lineages <- function(obs_genes_ls, meta_genes, avg_gene_exp, output_dir, datatag){
  # obs_genes_ls <- c('CDK14','HMGB2','TOP2A','PABPC1')
  # obs_genes_ls <- c('CDK1','MAEL','CDK14','FOXP1','CTSZ')
  # obs_genes_ls <- c('TCF4')
    
  # unique(avg_gene_exp$gene)
  for(i in seq(obs_genes_ls)){
    gsymb <- obs_genes_ls[i]
    obs_gene <- meta_genes[gsymb,'ens_gene_id']  #obs_gene: ens genes id
    
    if(!is.null(obs_gene)){
      sigf_gene_exp <- avg_gene_exp %>% 
        dplyr::filter(gene==obs_gene)
      print(dim(sigf_gene_exp))
      # p <- viz_smooth_genes_cluster(gsymb, sigf_gene_exp, 
      #                               obs_gene, paste0(gsymb, ': ',names(gsymb)), output_dir,
      #                               verbose_legend = 'none')
      verbose_lg <- c(0.3, 0.75)
      # verbose_lg
      p <- viz_smooth_genes_cluster_v2(datatag, gsymb, sigf_gene_exp, 
                                  obs_gene, paste0(gsymb), output_dir, #, ': ',names(gsymb)
                                  verbose_legend = verbose_lg)
      saveRDS(p, paste0(output_dir,gsymb,'.rds'))
      # p1 <- viz_smooth_genes_cluster(paste0(gsymb,'_lg'), sigf_gene_exp, 
      #                               obs_gene, paste0(gsymb, ': ',names(gsymb)), output_dir,
      #                               verbose_legend = 'right')
      p1 <- viz_smooth_genes_cluster_v2(datatag, gsymb, sigf_gene_exp, 
                                       obs_gene, paste0(gsymb), output_dir,
                                       verbose_legend = 'none')
      saveRDS(p1, paste0(output_dir,gsymb,'_without_lg.rds')) 
    }
  }
  
  
}


# We find genes at the top that are also ranked as DE for the differentiated cell type. 
# What is especially interesting are genes that have different expression patterns 
# but no different expression at the differentiated cell type level. 
# We therefore sort the genes according to the sum of square of their rank 
# in increasing Wald statistics for the patternTest 
# and their rank in decreasing/increasing Wald statistics for the diffEndTest or earlyDETest
# patternRes: matrix of gene and statistical values of patternTest
# suppTestRes: matrix of gene and statistical values of diffEndTest, or earlyDETest
get_square_rank <- function(patternRes, suppTestRes, save_dir, base_name='patternTest_diffEndTest',
                            xtitle="patternTest Wald Stat (log)", 
                            ytitle="diffEndTest Wald Stat (log)")
{
  
  # Get only significant genes from patternTest
  patternRes <- patternRes %>%
    dplyr::filter(pvalue<0.05)%>%
    dplyr::rename(pattern=waldStat, pvalue_pattern=pvalue)%>%
    dplyr::select(pattern, pvalue_pattern)%>%
    dplyr::arrange(desc(pattern))%>%
    tibble::rownames_to_column('ens_gene_id')
  print(head(patternRes))
  print(dim(patternRes))
  
  suppTestRes <- suppTestRes %>%
    dplyr::filter(pvalue<0.05)%>%
    dplyr::rename(statSecondTest=waldStat, pvalue_secondTest=pvalue)%>%
    dplyr::select(statSecondTest, pvalue_secondTest)%>%
    tibble::rownames_to_column('ens_gene_id')
  print(head(suppTestRes))
  print(dim(suppTestRes))
  
  compare_df <- patternRes %>% inner_join(suppTestRes, by=c('ens_gene_id'))
  print(dim(compare_df))
  compare_df$transientScore <- rank(compare_df$statSecondTest, ties.method = "random")^2 + rank(compare_df$pattern, ties.method = "random")^2
  
  
  # Add annotation info
  # ref_genes_fn <- '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/Symbol_ensembl.csv'
  # meta_genes <- data.table::fread(ref_genes_fn) %>% as.data.frame()
  # meta_genes <- meta_genes %>%
  #   dplyr::rename(ens_gene_id=Ensembl, gene_symbol=Symbol)
  # compare_df <- compare_df %>% inner_join(meta_genes, by=c('ens_gene_id'))
  annots <- annotables::grch38 %>% 
    # dplyr::select(gene_symbol = symbol, chr, start, end) %>% 
    dplyr::select(ens_gene_id = ensgene, gene_symbol=symbol, chr, description)
  annots <- annots[!duplicated(annots),]
  compare_df <- compare_df %>% left_join(annots, by='ens_gene_id')
  print(dim(compare_df))
  # patternRes <- patternRes %>% left_join(annots, by='ens_gene_id')
  compare_df <- compare_df[order(compare_df$transientScore, decreasing = T),]
  data.table::fwrite(compare_df, paste0(save_dir,base_name, ".csv"))
  print(dim(compare_df))
  # compare_df <- compare_df[order(compare_df$pattern, decreasing = T),]
  # print(head(compare_df[,1:8]))
  # p <- ggplot(compare_df, aes(x = log(pattern), y = log(statSecondTest))) +
  #   geom_point(aes(col = transientScore)) +
  #   labs(x = xtitle, y = ytitle) +
  #   scale_color_continuous(low = "blue", high = "red") +
  #   theme_classic()
  
  # png(paste0(save_dir,base_name, ".png"), height = 2*350, width=2*450, res = 2*72)
  # print(p)
  # dev.off()
  
  # res <- list(p=p, stat=compare_df)
  # saveRDS(res, paste0(save_dir,base_name, ".rds"))
  
  # return(res)
  return(compare_df)
}

# earlyDE genes between lineage or endDETest genes between lineages
get_significant_genes_drug_treatment <- function(patternRes, suppTestRes){
  # Get only significant genes from patternTest
  patternRes <- patternRes %>%
    dplyr::filter(pvalue<0.05)%>%
    dplyr::rename(pattern=waldStat, pvalue_pattern=pvalue)%>%
    dplyr::select(pattern, pvalue_pattern)%>%
    tibble::rownames_to_column('ens_gene_id')
  # print(head(patternRes))
  print(dim(patternRes))
  
  # suppTestRes <- suppTestRes %>%
  #   dplyr::filter(pvalue<0.05)%>%
  #   dplyr::rename(statSecondTest=waldStat, pvalue_secondTest=pvalue)%>%
  #   # dplyr::select(statSecondTest, pvalue_secondTest)%>%
  #   tibble::rownames_to_column('ens_gene_id')
  # # print(head(suppTestRes))
  # print(dim(suppTestRes))
  
  compare_df <- patternRes %>% inner_join(suppTestRes, by=c('ens_gene_id'))
  print(dim(compare_df))
  # Add annotation info
  annots <- annotables::grch38 %>% 
    # dplyr::select(gene_symbol = symbol, chr, start, end) %>% 
    dplyr::select(ens_gene_id = ensgene, gene_symbol=symbol, chr, description)
  annots <- annots[!duplicated(annots),]
  compare_df <- compare_df %>% left_join(annots, by='ens_gene_id')
  compare_df <- compare_df[order(compare_df$pattern, decreasing = T),]
  print(dim(compare_df))
  return(compare_df)
}


# verbose_legend: 'left','right', or x, y position ex: c(0.1, 0.8)
viz_smooth_genes_cluster_v2 <- function(datatag, gene_module, sigf_gene_exp, g, gsymb, 
                                        save_dir, verbose_legend='none'){
  
  res <- get_lineage_colors_linetype(sigf_gene_exp, datatag)
  sigf_gene_exp <- res$df
  unique(res$plt_setting$lineage)
  cols_use <- res$plt_setting$used_color
  names(cols_use) <- res$plt_setting$lineage
  linetypes <- res$plt_setting$linetypes
  names(linetypes) <- res$plt_setting$lineage
  line_sizes <- res$plt_setting$line_sizes
  names(line_sizes) <- res$plt_setting$lineage
  
  my_font <- "Helvetica"
  df <- sigf_gene_exp[sigf_gene_exp$gene == sigf_gene_exp$gene[1],]
  
  p <- ggplot(data = df, aes(x = time, y = yhat)) +
    geom_point(alpha = 0) +
    labs(title = gsymb,  x = "Pseudotime", y = "Relative Expression") + #paste0("Cluster ", xx)
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, family=my_font, size=12),
          axis.text = element_text(size=8, hjust = 0.5, family=my_font),
          legend.text = element_text(size=11, family=my_font),
          legend.position = verbose_legend,
          legend.title = element_blank(),
          legend.background=element_rect(fill = alpha("white", 0)))
  if('lineage_desc' %in% colnames(sigf_gene_exp)){
    sigf_gene_exp$lineage_desc <- gsub('Lineage ','L',sigf_gene_exp$lineage_desc)  
  }else{
    sigf_gene_exp$lineage_desc <- sigf_gene_exp$lineage
  }
  
  # unique(sigf_gene_exp$lineage_desc)
  for(g in unique(sigf_gene_exp$gene)){
    gene_exp <- sigf_gene_exp[sigf_gene_exp$gene == g,]
    p <- p + geom_line(data = gene_exp,
                       aes(x=time, y=yhat, col = lineage_desc, group = lineage_desc, size=lineage_desc), lwd = 1) + 
      scale_color_manual(values = cols_use)+
      # scale_linetype_manual(values=linetypes)+
      scale_size_manual(values=line_sizes)
  }
  
  # p <- p + scale_color_manual(values = c("orange", "darkseagreen3",'blue'))  
  png(paste0(save_dir, "gene_module_",gene_module, ".png"), height = 2*250, width=2*350, res = 2*72)
  print(p)
  dev.off()
  return(p)
}  

viz_smooth_genes_cluster <- function(gene_module, sigf_gene_exp, g, gsymb, save_dir, verbose_legend='none'){
  my_font <- "Helvetica"
  # t <- sigf_gene_exp %>%
  #   dplyr::filter(module==gene_module)
  t <- sigf_gene_exp
  # print(dim(t))
  df <- sigf_gene_exp[sigf_gene_exp$gene == sigf_gene_exp$gene[1],]
  
  p <- ggplot(data = df, aes(x = time, y = yhat)) +
    geom_point(alpha = 0) +
    labs(title = gsymb,  x = "Pseudotime", y = "Log Normalized Expr") + #paste0("Cluster ", xx)
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, family=my_font, size=12),
          axis.text = element_text(size=8, hjust = 0.5, family=my_font),
          legend.position = verbose_legend)
  
  for(g in unique(t$gene)){
    gene_exp <- sigf_gene_exp[sigf_gene_exp$gene == g,]
    p <- p + geom_line(data = gene_exp,
                       aes(x=time, y=yhat, col = lineage_desc, group = lineage_desc), lwd = 0.5)
  }
  # p <- p + scale_color_manual(values = c("orange", "darkseagreen3",'blue'))  
  png(paste0(save_dir, "gene_module_",gene_module, ".png"), height = 2*250, width=2*350, res = 2*72)
  print(p)
  dev.off()
  return(p)
}  
viz_smooth_gene_exp <- function(sigf_gene_exp, g, gsymb, save_dir){
  gene_exp <- sigf_gene_exp[sigf_gene_exp$gene == g,]
  # dim(gene_exp)
  p <- ggplot(data = gene_exp, aes(x = time, y = yhat)) +
    geom_line(aes(col = lineage_desc, group = lineage_desc), lwd = 1.5) + 
    labs(title = paste0("Gene ", gsymb),  x = "Pseudotime", y = "Normalized expression") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  
  png(paste0(save_dir, gsymb, ".png"), height = 2*350, width=2*450, res = 2*72)
  print(p)
  dev.off()
}
get_signficant_gene_reversibility_SA535 <- function(output_dir){
  datatag <- 'SA535'
  input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')
  output_dir <- paste0(save_dir,'tradeSeq_3000/')
  pvalue_thrs <- 0.05
  # only for SA535, end test
  endRes <- readRDS(paste0(output_dir, "endRes_out.rds"))
  patternRes <- readRDS(paste0(output_dir, "patternRes_out.rds"))
  earlyDERes <- readRDS(paste0(output_dir, "earlyDERes_out.rds"))
  ts_sce <- readRDS(paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/tradeSeq_3000/', "fitGAM_out.rds"))
  
  # View(head(patternRes))
  # Get reversibility genes
  suppTestRes <- endRes %>%
    dplyr::select(waldStat_1vs2, pvalue_1vs2) %>%
    dplyr::rename(waldStat=waldStat_1vs2, pvalue=pvalue_1vs2)
  
  suppTestRes <- earlyDERes %>%
    dplyr::select(waldStat_1vs2, pvalue_1vs2) %>%
    dplyr::rename(waldStat=waldStat_1vs2, pvalue=pvalue_1vs2)
  dim(earlyDERes)
  dim(suppTestRes)
  head(suppTestRes)
  summary(suppTestRes$waldStat)
  ytitle="diffEndTest Wald Stat (log)"
  xtitle="patternTest Wald Stat (log)"
  base_name='patternTest_diffEndTest_SA535_l2_Rx_lg1_RxH'
  base_name='patternTest_earlyDERes_SA535_l2_Rx_lg1_RxH'
  save_dir <- paste0(output_dir,'tradeSeq_SA535/reversibility_genes/')
  dir.create(save_dir)
  compare_df <- get_square_rank(patternRes, suppTestRes, save_dir, base_name, xtitle, ytitle)
  dim(compare_df)

  # summary(compare_df$statSecondTest)
  # sum(compare_df$statSecondTest>=40)
  # Get all significant genes
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/')
  save_dir <- paste0(save_dir,'tradeSeq_3000/tradeSeq_SA535/Rx_vs_UnRx/')
  dir.create(save_dir)
  dim(endRes)
  colnames(endRes)
  
  # suppTestRes <- endRes %>%
  #   dplyr::filter(pvalue_1vs4<pvalue_thrs | pvalue_2vs4<pvalue_thrs | pvalue_3vs4<pvalue_thrs | pvalue_1vs2 < pvalue_thrs)# 
  suppTestRes <- endRes %>%
    dplyr::filter(pvalue_2vs4<pvalue_thrs)# 
  dim(suppTestRes)
  base_name='SA535_Rx_vs_UnRx_endTest'
  suppTestRes <- suppTestRes %>%
    dplyr::rename(waldStat_global_endTest=waldStat, pvalue_global_endTest=pvalue)%>%
    dplyr::select(-df)%>%
    tibble::rownames_to_column('ens_gene_id')
  dim(patternRes)
  stat_df <- get_significant_genes_drug_treatment(patternRes, suppTestRes)
  data.table::fwrite(stat_df, paste0(save_dir, base_name, '.csv'))
  dim(stat_df)
  
  # lineage 1, 2, 3: Rx, RxH, lineage 4: UnRx
  suppTestRes <- earlyDERes %>%
    dplyr::filter(pvalue_2vs4<pvalue_thrs)#
  base_name='SA535_Rx_vs_UnRx_earlyTest'
  dim(suppTestRes)
  colnames(earlyDERes)
  colnames(patternRes)
  suppTestRes <- suppTestRes %>%
    dplyr::rename(waldStat_global_earlyTest=waldStat, pvalue_global_earlyTest=pvalue)%>%
    dplyr::select(-df)%>%
    tibble::rownames_to_column('ens_gene_id')
  stat_df <- get_significant_genes_drug_treatment(patternRes, suppTestRes)
  data.table::fwrite(stat_df, paste0(save_dir, base_name, '.csv'))
  # View(head(stat_df))
  
  # Visualization
  # yhatSmooth <- predictSmooth(ts_sce, gene = mockGenes, nPoints = 50, tidy = FALSE)
  # heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))),
  #                        cluster_cols = FALSE,
  #                        show_rownames = FALSE, 
  #                        show_colnames = FALSE)
  # 
  # 
  # obs_gene <- compare_df$ens_gene_id[1:20]
  ysmooth <-  tradeSeq::predictSmooth(models = ts_sce,
                                      gene = obs_genes,
                                      nPoints = 100)
  # dim(ysmooth)
  # 
  # class(ysmooth)
  # # ysmooth <- as.data.frame(ysmooth)
  # unique(ysmooth$lineage)
  
  # preprocess_mat <- as.matrix(counts(ts_sce[obs_gene,]))
  # gene_module_df <- get_genes_modules(preprocess_mat=preprocess_mat, 
  #                                     resolution=0.5)
  # data.table::fwrite(gene_module_df, paste0(output_dir,'gene_modules_',base_name,'.csv'))
  # dim(gene_module_df)
  # summary(as.factor(gene_module_df$module))
  # unique(gene_module_df$supermodule)
  # colnames(gene_module_df)
  # gene_module_df$id[1]
  # colnames(ysmooth)
  # ysmooth$gene[1]
  # unique(sigf_gene_exp$time)
  # sigf_gene_exp <- ysmooth %>% inner_join(gene_module_df, by=c('gene'='id'))
  # dim(sigf_gene_exp)
  sigf_gene_exp <- ysmooth
  # summary(as.factor(sigf_gene_exp$module))
  sigf_gene_exp$lineage <- paste0('Lineage ',sigf_gene_exp$lineage)
  
  sigf_gene_exp <- sigf_gene_exp %>%
    # dplyr::filter(lineage %in% compare_df$ens_gene_id[1:50])
    dplyr::filter(lineage %in% c("Lineage 1","Lineage 2","Lineage 4"))
  
  unique(sigf_gene_exp$lineage)
  meta_lineage <- data.frame(lineage=c("Lineage 1","Lineage 2","Lineage 4"),
                             lineage_desc=c("L1: RxH,Rx","L2: Rx, RxH","L4: UnRx"))
  sigf_gene_exp <- sigf_gene_exp %>% inner_join(meta_lineage, by=c('lineage'))
  # gene_module <- 1
  # viz_smooth_genes_cluster(gene_module, sigf_gene_exp, g, gsymb, save_dir)
  
  # Visualize given gene
  datatag <- 'SA535'
  save_dir <- output_dir
  # obs_gene_ls <- c('MYC','FOS')
  # obs_gene_ls <- c('MYC','HIST1H4C','MKI67','ODAM')
  
  obs_gene_ls <- c('ODAM')
  obs_gene_ls <- c('CD99','ODAM')
  obs_gene_ls <- c('S100A6','SAA1')
  for(gsymb in obs_gene_ls){
    # grep('CK',compare_df$gene_symbol,value=T)
    # t <- compare_df %>%
    #   dplyr::filter(gene_symbol==gsymb)
    t <- meta_genes %>%
      dplyr::filter(gene_symbol == gsymb)
    if(dim(t)[1]>0){
      g <- t$ens_gene_id
      viz_smooth_gene_exp(sigf_gene_exp, g, paste0(datatag, ": ",gsymb), save_dir)
    }
    
  }
  
}  
get_signficant_gene_reversibility_SA609 <- function(output_dir){
  # only for SA609
  pvalue_thrs <- 0.05
  datatag <- 'SA609'
  base_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/'
  input_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/')
  save_dir <- paste0('/home/htran/storage/datasets/drug_resistance/rna_results/',datatag,'_rna/slingshot_trajectory/withBE_SA609_v2/')
  nfeatures_use <- 3000
  output_dir <- paste0(save_dir,'tradeSeq_3000/')
  
  earlyDERes <- readRDS(paste0(output_dir, "earlyDERes_out.rds"))
  patternRes <- readRDS(paste0(output_dir, "patternRes_out.rds"))
  endRes <- readRDS(paste0(output_dir, "endRes_out.rds"))
  
  # get signif gene between Rx, RxH vs. UnRx
  genes_use <- endRes %>%
    dplyr::filter(pvalue_1vs3 < pvalue_thrs | pvalue_2vs3 < pvalue_thrs) %>%
    tibble::rownames_to_column('ens_gene_id') %>%
    dplyr::pull(ens_gene_id)
  # genes_use <- earlyDERes %>%
  #   dplyr::filter(pvalue_1vs3<0.05 | pvalue_2vs3 < 0.05) %>%
  #   tibble::rownames_to_column('ens_gene_id') %>%
  #   dplyr::pull(ens_gene_id)
  print(length(genes_use))
  # earlyDERes <- earlyDERes[genes_use,]
  endRes <- endRes[genes_use,]
  
  # suppTestRes <- earlyDERes %>%
  #   dplyr::select(waldStat_1vs2, pvalue_1vs2) %>%
  #   dplyr::rename(waldStat=waldStat_1vs2, pvalue=pvalue_1vs2)
  suppTestRes <- endRes %>%
    dplyr::select(waldStat_1vs2, pvalue_1vs2) %>%
    dplyr::rename(waldStat=waldStat_1vs2, pvalue=pvalue_1vs2)
  dim(suppTestRes)
  head(suppTestRes)
  # ytitle="earlyDETest Wald Stat (log)"
  ytitle="diffEndTest Wald Stat (log)"
  xtitle="patternTest Wald Stat (log)"
  base_name='patternTest_diffEndTest_SA609_l1_Rx_lg2_RxH'
  save_dir <- paste0(output_dir,'SA609_tradeSeq/reversibility_genes/')
  res <- get_square_rank(patternRes, suppTestRes, save_dir, base_name, xtitle, ytitle)
  compare_df <- res$stat
  
  
  # Get all significant genes
  save_dir <- paste0(output_dir,'SA609_tradeSeq/Rx_vs_UnRx/')
  dim(endRes)
  colnames(endRes)
  genes_use <- endRes %>%
    dplyr::filter(pvalue_1vs3 < pvalue_thrs | pvalue_2vs3 < pvalue_thrs) %>%
    tibble::rownames_to_column('ens_gene_id') %>%
    dplyr::pull(ens_gene_id)
  print(length(genes_use))
  endRes <- endRes[genes_use,]
  suppTestRes <- endRes
  base_name='SA609_Rx_vs_UnRx_endTest'
  dim(suppTestRes)
  colnames(suppTestRes)
  colnames(patternRes)
  suppTestRes <- suppTestRes %>%
    dplyr::rename(waldStat_global_endTest=waldStat, pvalue_global_endTest=pvalue)%>%
    dplyr::select(-df)%>%
    tibble::rownames_to_column('ens_gene_id')
  stat_df <- get_significant_genes_drug_treatment(patternRes, suppTestRes)
  data.table::fwrite(stat_df, paste0(save_dir, base_name, '.csv'))
  
  # Get all significant genes
  genes_use <- earlyDERes %>%
    dplyr::filter(pvalue_1vs3 < pvalue_thrs | pvalue_2vs3 < pvalue_thrs) %>%
    tibble::rownames_to_column('ens_gene_id') %>%
    dplyr::pull(ens_gene_id)
  print(length(genes_use))
  earlyDERes <- earlyDERes[genes_use,]
  suppTestRes <- earlyDERes
  base_name='SA609_Rx_vs_UnRx_earlyTest'
  dim(suppTestRes)
  colnames(suppTestRes)
  colnames(patternRes)
  suppTestRes <- suppTestRes %>%
    dplyr::rename(waldStat_global_earlyTest=waldStat, pvalue_global_earlyTest=pvalue)%>%
    dplyr::select(-df)%>%
    tibble::rownames_to_column('ens_gene_id')
  stat_df <- get_significant_genes_drug_treatment(patternRes, suppTestRes)
  data.table::fwrite(stat_df, paste0(save_dir, base_name, '.csv'))
  # View(head(stat_df))
  
  
  # Visualization
  obs_gene <- compare_df$ens_gene_id[1:50]
  ysmooth <-  tradeSeq::predictSmooth(models = ts_sce, 
                                      gene = obs_gene, 
                                      nPoints = 100)
  dim(ysmooth)
  
  class(ysmooth)
  # ysmooth <- as.data.frame(ysmooth)
  unique(ysmooth$lineage)
  
  preprocess_mat <- as.matrix(counts(ts_sce[obs_gene,]))
  gene_module_df <- get_genes_modules(preprocess_mat=preprocess_mat, 
                                      resolution=0.5)
  data.table::fwrite(gene_module_df, paste0(output_dir,'gene_modules_',base_name,'.csv'))
  dim(gene_module_df)
  summary(as.factor(gene_module_df$module))
  unique(gene_module_df$supermodule)
  colnames(gene_module_df)
  gene_module_df$id[1]
  colnames(ysmooth)
  ysmooth$gene[1]
  unique(sigf_gene_exp$time)
  sigf_gene_exp <- ysmooth %>% inner_join(gene_module_df, by=c('gene'='id'))
  dim(sigf_gene_exp)
  summary(as.factor(sigf_gene_exp$module))
  # sigf_gene_exp$lineage <- factor(sigf_gene_exp, levels=sort(unique(sigf_gene_exp$lineage)))
  sigf_gene_exp$lineage <- paste0('Lineage ',sigf_gene_exp$lineage)
  
  # sigf_gene_exp <- sigf_gene_exp %>%
  #   # dplyr::filter(lineage %in% compare_df$ens_gene_id[1:50])
  #   dplyr::filter(lineage %in% c("Lineage 1","Lineage 2"))
  
  
  
  unique(sigf_gene_exp$lineage)
  meta_lineage <- data.frame(lineage=c("Lineage 1","Lineage 2","Lineage 3"),
                             lineage_desc=c("Lineage 1: Rx","Lineage 2: RxH","Lineage 3: UnRx"))
  sigf_gene_exp <- sigf_gene_exp %>% inner_join(meta_lineage, by=c('lineage'))
  gene_module <- 1
  viz_smooth_genes_cluster(gene_module, sigf_gene_exp, g, gsymb, save_dir)
  
  # Visualize given gene
  g <- obs_gene[1]
  gsymb <- compare_df$gene_symbol[1]
  # grep('CK',compare_df$gene_symbol,value=T)
  # gsymb <- 'MKI67'
  # gsymb <- 'EGFR'
  # gsymb <- 'CKS1B'
  # t <- compare_df %>%
  #   dplyr::filter(gene_symbol==gsymb)
  # dim(t)
  # g <- t$ens_gene_id
  viz_smooth_gene_exp(sigf_gene_exp, g, gsymb, save_dir)
  
}

# Not a good algo but ok 
convert_range_pseudotime <- function(x, ref){
  # x <- c(1,2,5)
  # ref <- c(1.5,2.2, 5.1, 6)
  mx <- max(ref)
  pseudo_2 <- c()
  for(i in x){
    min_dis <- mx
    newval <- -1
    for(r in ref){
      if(abs(r-i)<min_dis){
        min_dis <- abs(r-i)
        newval <- r
      }
    }
    pseudo_2 <- c(pseudo_2,newval)
  }
  return(pseudo_2)
}

# ts_sce: tradeSeq obj
# obs_genes_df: genes_use, gene_symbol
# meta_genes <- data.frame(genes_use=patternRes$ens_gene_id,
#                          gene_symbol=patternRes$gene_symbol, row.names = patternRes$gene_symbol)
viz_heatmap_transient_genes_SA609 <- function(ts_sce, obs_genes_df, meta_genes){
  dh_lineage <- 2
  save_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/"
  datatag <- 'SA609'
  sce <- readRDS(paste0(save_dir, datatag,'_3000_rd_sce.rds'))
  meta_cells <- colData(sce)
  unique(meta_cells$clone)
  obs_clones_dh <- c('B','C','C_D','D')
  dh_cells <- meta_cells %>%
    as.data.frame()%>%
    dplyr::filter(treatmentSt=='UTU' & clone %in% obs_clones_dh)%>%
    dplyr::pull(cell_id)
  length(dh_cells)
  
  # get tradeSeq info corresponding to drug holiday lineage
  dm <- colData(ts_sce)$tradeSeq$dm # design matrix
  X <- colData(ts_sce)$tradeSeq$X # linear predictor
  slingshotColData <- colData(ts_sce)$slingshot
  
  pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
                                       pattern = "pseudotime")]
  if (is.null(dim(pseudotime))) {
    pseudotime <- matrix(pseudotime, ncol = 1) 
  }
  dh_lg <- paste0('pseudotime.Lineage',dh_lineage)
  pseudotime <- pseudotime %>%
    as.data.frame()%>%
    dplyr::select(!!sym(dh_lg))
  pseudotime <- pseudotime[dh_cells,, drop=F]
  
  # pseudotime$pseudo_2 <- as.numeric(ggplot2::cut_number(unique(pseudotime$pseudotime.Lineage2)),100)
  
  pseudotime$pseudotime.Lineage2[1:100]
  dim(pseudotime)
  summary(pseudotime$pseudotime.Lineage2)
  yhatSmooth$lineage[1]
  t <- yhatSmooth %>%
    as.data.frame()%>%
    dplyr::filter(lineage==dh_lineage)
  dim(t)
  t$time[1:100]
  length(unique(pseudotime$dh_time))
  pseudotime$dh_time <- convert_range_pseudotime(pseudotime$pseudotime.Lineage2, unique(t$time))
  pseudotime$dh_time
  
  pseudo <- unique(t$time)
  labels <- ifelse(pseudo %in% unique(pseudotime$dh_time),'RxH','No_RxH')
  anno_col <- ComplexHeatmap::HeatmapAnnotation(treatment=factor(c(rep('No_RxH',nb_estimated_points),
                                                                   labels,
                                                                   rep('No_RxH',nb_estimated_points))),
                                                col = list(treatment = c(No_RxH = "#D0D0D0", RxH = "#0D4848")))
  
  nb_genes <- dim(obs_genes_df)[1]
  nb_estimated_points <- 100
  yhatSmooth <- tradeSeq::predictSmooth(ts_sce, gene = obs_genes_df$genes_use, nPoints = nb_estimated_points)
  yhatSmooth$desc <- paste0(yhatSmooth$lineage,'_', yhatSmooth$time)
  l1 <- yhatSmooth %>%
    # dplyr::arrange(lineage)%>%
    dplyr::select(-lineage, -time)%>%
    tidyr::pivot_wider(names_from = 'desc', values_from = 'yhat')%>%
    tibble::column_to_rownames('gene')
  obs_cells_df <- data.frame(lineage=c(rep('L1',nb_estimated_points),rep('L2',nb_estimated_points),rep('L3',nb_estimated_points)))
  rownames(obs_cells_df) <- colnames(l1)
  nb_estimated_points <- 100
  lgs <- unique(yhatSmooth$lineage)
  lineage <- c()
  for(lg in lgs){
    lineages <- c(lineage, rep(paste0('Lineage',lg), nb_estimated_points))
  }
  obs_cells_df <- data.frame(lineage=lineages)
  
  max_exp <- 3
  tmp <- t(scale(t(l1)))
  tmp[tmp>max_exp] <- max_exp
  tmp[tmp<-max_exp] <- -max_exp
  rownames(obs_cells_df) <- colnames(l1)
  rownames(tmp) <- meta_genes[rownames(tmp),'gene_symbol']
  srn <- ifelse(nb_genes<=50,T,F)
  # anno_col <- NULL
  # anno_col <- ComplexHeatmap::HeatmapAnnotation(treatment=factor(lineages))
  p <- ComplexHeatmap::Heatmap(tmp, na_col = "white",
                               show_column_names=F,
                               show_row_names = srn,
                               cluster_rows=T,cluster_columns=F,
                               name = "Avg Exp", 
                               # row_order = sort(rownames(test)),
                               row_split= obs_genes_df,
                               row_title_rot = 0,
                               row_gap = unit(2, "mm"),
                               column_split = obs_cells_df, 
                               # column_title = "ddd",
                               column_gap = unit(2, "mm"),
                               column_names_gp = grid::gpar(fontsize = 6),
                               row_names_gp = grid::gpar(fontsize = 7),
                               show_heatmap_legend = T,
                               top_annotation=anno_col,
                               # left_annotation = left_anno,
                               # cell_fun = cell_func,
                               row_dend_reorder=F,
                               column_title = "Pseudotime", 
                               column_title_side = "bottom")
  # add title
  hz <- ifelse(nb_genes<=50,(nb_genes*8+10),550)
  png(paste0(output_dir,"clusters_genes.png"), height = 2*hz, width=2*750, res = 2*72)
  print(p)
  dev.off()
  return(p)
}

write_csv_genes_modules <- function(obs_genes, datatag, tag, output_dir){
  
  annots <- annotables::grch38 %>% 
    # dplyr::select(gene_symbol = symbol, chr, start, end) %>% 
    dplyr::select(ens_gene_id = ensgene, gene_symbol=symbol, chr, description)
  annots <- annots[!duplicated(annots$ens_gene_id),]
  # meta_genes <- data.frame(ens_gene_id=gene_module_df$id)
  obs_genes_df <- data.frame(ens_gene_id=obs_genes)
  obs_genes_df <- obs_genes_df %>% left_join(annots, by=c('ens_gene_id'))
  data.table::fwrite(obs_genes_df, paste0(output_dir, datatag,'_', tag, '.csv'))
}  
get_meta_genes <- function(obs_genes_df){
  # obs_genes_df <- total_genes
  # meta_genes <- data.frame(ens_gene_id=gene_module_df$id)
  # meta_genes <- data.frame(ens_gene_id=obs_genes_df$genes_use)
  obs_genes_df <- as.data.frame(obs_genes_df)
  if('ens_gene_id' %in% colnames(obs_genes_df)){
    meta_genes <- obs_genes_df
  }else{
    meta_genes <- obs_genes_df %>%
      dplyr::rename(ens_gene_id=genes_use)
  }
  if(!'gene_symbol' %in% colnames(obs_genes_df)){
    annots <- annotables::grch38 %>% 
      # dplyr::select(gene_symbol = symbol, chr, start, end) %>% 
      dplyr::select(ens_gene_id = ensgene, gene_symbol=symbol, chr, description)
    annots <- annots[!duplicated(annots$ens_gene_id),]
    meta_genes <- meta_genes %>% 
      # as.data.frame %>% 
      dplyr::left_join(annots, by=c('ens_gene_id'))  
  }
  rownames(meta_genes) <- meta_genes$ens_gene_id
  # meta_genes$gene_symbol[1]
  # dim(meta_genes)
  # meta_genes <- data.frame(genes_use=patternRes$ens_gene_id,
  #                          gene_symbol=patternRes$gene_symbol, row.names = patternRes$gene_symbol)
  # head(meta_genes)
  return(meta_genes)
}

run_tradeSeq <- function(counts, pseudotime, cellWeights, output_dir){
  ts_sce <- tradeSeq::fitGAM(counts = counts,
                             pseudotime = pseudotime, cellWeights = cellWeights,
                             nknots = 12, verbose = FALSE,
                             parallel=T,
                             BPPARAM = BiocParallel::MulticoreParam(workers = 5))
  saveRDS(ts_sce, paste0(output_dir, "fitGAM_out.rds"))
  # assoRes <- associationTest(ts_sce)
  # saveRDS(assoRes, paste0(output_dir, "assoRes_out.rds"))
  # assoRes <- readRDS(paste0(output_dir, "assoRes_out.rds"))
  
  startRes <- startVsEndTest(ts_sce, lineages=TRUE)
  saveRDS(startRes, paste0(output_dir, "startRes_out.rds"))
  # print(head(startRes))
  
  endRes <- diffEndTest(ts_sce, pairwise=TRUE)
  # print(head(endRes))
  saveRDS(endRes, paste0(output_dir, "endRes_out.rds"))
  
  patternRes <- patternTest(ts_sce)
  # oPat <- order(patternRes$waldStat, decreasing = TRUE)
  # print(head(rownames(patternRes)[oPat]))
  saveRDS(patternRes, paste0(output_dir, "patternRes_out.rds"))
  
  
  earlyDERes <- tradeSeq::earlyDETest(ts_sce, l2fc = 0.5, pairwise=TRUE)
  saveRDS(earlyDERes, paste0(output_dir, "earlyDERes_out.rds"))
  dim(earlyDERes)
  
}
