suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot, quietly = TRUE)
  # library(mclust, quietly = TRUE)
  library(grDevices)
  library(RColorBrewer)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(monocle3)
  library(circlize)
  library(ComplexHeatmap)
})
suppressPackageStartupMessages(library(circlize))
options(dplyr.summarise.inform = FALSE)
options(tidyverse.quiet = TRUE)
library(extrafont)
# font_import(prompt=F, paths ='/usr/share/fonts/truetype/myfonts/') # import Helvetica font
fonts()
my_font <- "Helvetica"


get_slingshot_pseudotime_v2 <- function(sce, save_dir, datatag, 
                                        start_cls=NULL, cl=NULL, rd_use='PCA'){
  save_dir <- paste0(save_dir,'slingshot_output/')
  dir.create(save_dir, showWarnings=F)
  rd_use='PCA'
  # UU in case of SA535, and U in other case, starting points - seeding cluster
  meta_info <- SingleCellExperiment::colData(sce) %>% as.data.frame()
  start_cls <- get_start_cluster(meta_info, datatag, root_ts=c('U'))
  print(start_cls)
  rd <- SingleCellExperiment::reducedDims(sce)[[rd_use]]
  umaps <- SingleCellExperiment::reducedDims(sce)[['UMAP']]
  dim(umaps)
  print(dim(rd))
  class(rd)
  if(is.null(cl)){
    cl <- colData(sce)[,'cluster_label'] # Gaussian mixture model clustering - optimize nb clusters  
  }
  # print(unique(cl))
  # length(cl)  
  lin1 <- slingshot::getLineages(rd, cl, start.clus = as.numeric(start_cls))
  approx_points_use <- 400
  crv1 <- slingshot::getCurves(lin1, approx_points = approx_points_use)
  class(crv1)
  saveRDS(crv1, paste0(save_dir, "slingshot_pseudotime_",datatag,'_',paste(start_cls, collapse='_'),'_',rd_use,"_crv.rds"))
  
  pseudotime <- slingPseudotime(crv1, na = FALSE) %>% as.data.frame()
  dim(pseudotime)
  pseudotime$cell_id <- rownames(pseudotime)
  # View(head(pseudotime))
  cellWeights <- slingshot::slingCurveWeights(crv1) %>% as.data.frame()
  dim(cellWeights)
  cellWeights$cell_id <- rownames(cellWeights)
  data.table::fwrite(pseudotime, paste0(save_dir, "slingshot_",datatag,'_',paste(start_cls, collapse='_'),'_',rd_use,"_pseudotime.csv"))
  data.table::fwrite(cellWeights, paste0(save_dir, "slingshot_",datatag,'_',paste(start_cls, collapse='_'),'_',rd_use,"_cellWeights.csv"))

  
  # For visualization
  # View(head(umaps))
  sum(rownames(umaps)==colnames(sce))
  crv_umap_embed <- embedCurves(crv1, umaps)
  saveRDS(crv_umap_embed, paste0(save_dir, "slingshot_",datatag,'_',paste(start_cls, collapse='_'),"_UMAP_embed_crv.rds"))
  
  
  # res <- set_color_clusters(paste0('Cls_',cl))
  # lg <- plot_colors(res$clone_palette, save_dir, legend_label='clusters')
  # png(paste0(save_dir,"slingshot_",datatag,'_',paste(start_cls, collapse='_'),'_',rd_use,"_whole_lineages.png"), height = 2*500, width=2*600,res = 2*72)
  # plot(umaps, col = res$cols_use, asp = 1, pch = 16, cex = 1)
  # # lines(SlingshotDataSet(pt), type = 'lineages', lwd = 3, col = 'black')
  # lines(SlingshotDataSet(crv_umap_embed), type = 'lineages', lwd = 2.5, col = 'black')
  # dev.off()
  # 
  # if(datatag=='SA535'){
  #   sce$treatmentSt <- sce$treat  
  # }
  # 
  # # res2 <- set_color_treatmentSt(colData(sce)[,'treatmentSt'], datatag, cols=NULL)
  # res2 <- set_color_treatmentSt(colData(sce)[,'treatmentSt'], datatag, cols=NULL)
  # lg <- plot_colors(res2$clone_palette, save_dir, legend_label='treatmentSt')
  # png(paste0(save_dir,"slingshot_",datatag,'_',paste(start_cls, collapse='_'),'_',rd_use,"_lineage_treatmentSt.png"), height = 2*500, width=2*600,res = 2*72)
  # plot(umaps, col = res2$cols_use, asp = 1, pch = 16, cex = 1)
  # # lines(SlingshotDataSet(pt), type = 'lineages', lwd = 3, col = 'black')
  # lines(SlingshotDataSet(crv_umap_embed), type = 'lineages', lwd = 2.5, col = 'black')
  # dev.off()
  # 
  # plot_lineages(crv_umap_embed, umaps, save_dir)
  # plot_lineages_colorby_treatmentSt(crv_umap_embed, umaps, meta_info, datatag, save_dir)
  
} 

plot_whole_lineages_clones <- function(pt, rd, 
                                            curves, df_cls_centers, 
                                            output_dir, datatag, 
                                            downsample_ratio=0.25, cols_use=NULL,
                                            x_plt='UMAP_1', y_plt='UMAP_2'){
  col_plt <- 'clone'
  
  # pt <- pseudo_out 
  # rd_backup <- rd 
  # curves, df_cls_centers, 
  # save_figs_dir, datatag, cols_use=NULL,
  # x_plt='UMAP_1'
  # y_plt='UMAP_2'
  # head(curves)
  if(is.null(cols_use)){
    # res1 <- set_color_clones(rd$clone)
    # cols_use <- res1$clone_palette
    cols_use <- get_color_clones_v2(datatag, unique(rd$clone))
  }
  cols_use <- c(cols_use, 'NA'='#E8E8E8')
  
  ls_startp <- list()
  ls_endp <- list()
  for(lid in unique(curves$ligneage_id)){
    curves_tmp <- curves %>% 
      dplyr::filter(ligneage_id==lid)
    # print(dim(curves_tmp))
    obs_cls <- unique(c(curves_tmp$cls1, curves_tmp$cls2))
    df_cls_centers_tmp <- df_cls_centers %>% 
      dplyr::filter(cluster_label %in% obs_cls)
    # class(pt) dim(pt) colnames(pt)
    # print(dim(df_cls_centers_tmp))
    startp <- df_cls_centers_tmp %>% 
      dplyr::filter(cluster_label == curves_tmp$cls1[1])
    ls_startp[[curves_tmp$cls1[1]]] <- startp
    endp <- df_cls_centers_tmp %>% 
      dplyr::filter(cluster_label == curves_tmp$cls2[dim(curves_tmp)[1]])
    ls_endp[[curves_tmp$cls2[dim(curves_tmp)[1]]]] <- endp
    if(!lid %in% colnames(pt)){
      stop('Do not exist any results corresponding to given lineage, double check input data')
    }
  }
  
  # head(df_cls_centers)
  
  
  # View(head(pt))
  # head(rd)
  
  pt$cell_id <- rownames(pt)
  if(!'cell_id' %in% colnames(rd)){
    rd$cell_id <- rownames(rd)
  }
  if(downsample_ratio>0){
    print('Before downsampling')
    print(dim(rd))
    ## Original version
    # cells_use <- rd %>%
    #   # dplyr::filter(cluster_label %in% obs_cls)%>%
    #   dplyr::pull(cell_id)
    
    ## Downsample data, make plot clearer
    cells_use <- downsample_by_cluster_label_treatmentSt(rd, downsample_ratio=downsample_ratio) #0.15
    # rd$cell_id <- rownames(rd)
    # cells_use <- rd$cell_id
    print(length(cells_use))
    pt <- pt %>%
      dplyr::filter(cell_id %in% cells_use)
  }
  
  # rd <- rd %>% left_join(pt, by=c('cell_id'))
  rd <- rd %>% inner_join(pt, by=c('cell_id'))
  # print('After downsampling')
  print(dim(rd))
  # colors <- pal[cut(rd[,lid], breaks = 100)]
  # rd$treatmentSt <- ifelse(!is.na(rd[,lid]),rd$treatmentSt,'NA')
  
  # length(colors)
  df_cls_centers1 <- df_cls_centers # clusters center
  # Just to avoid overlap with cluster centers
  # df_cls_centers$x1 <- df_cls_centers$x1 + runif(dim(df_cls_centers)[1], -0.5, 0.5)
  # df_cls_centers$y1 <- df_cls_centers$y1 + runif(dim(df_cls_centers)[1], -0.5, 0.5)
  
  # Just to avoid overlap with cluster centers
  df_cls_centers$x1 <- df_cls_centers$x1 + 0.5
  df_cls_centers$y1 <- df_cls_centers$y1 - 0.2
  
  # Just for easy observe output in visualization
  if(datatag=='SA609'){
    # unique(df_cls_centers$cluster_label)
    meta_cluster <- data.frame(cluster_label=c('10','6','1','7','9','4','8','0','5','3','2'),
                               cls_lb=c('0','1','2','3','4','5','6','7','8','9','10'))
    df_cls_centers <- df_cls_centers %>% inner_join(meta_cluster, by=c('cluster_label'))
  }
  else if(datatag=='SA535'){
    meta_cluster <- data.frame(cluster_label=c('0','1','2','3','4','5','6','7','8','9','10'),
                               cls_lb=c('8','5','2','7','9','6','10','0','4','3','1'))
    df_cls_centers <- df_cls_centers %>% inner_join(meta_cluster, by=c('cluster_label'))
  }else if(datatag=='SA1035'){
    # meta_cluster <- data.frame(cluster_label=c('0','1','2','3','4','5','6','7','8'),
    #                            cls_lb=c('1','0','6','5','3','8','4','7','2')) # version 1 
    meta_cluster <- data.frame(cluster_label=c('0','1','10','6','3','7','9','5','11','4','2','8'),
                               cls_lb=c('0','1','2','3','4','5','6','7','8','9','10','11'))
    df_cls_centers <- df_cls_centers %>% inner_join(meta_cluster, by=c('cluster_label'))
    # df_cls_centers$cls_lb <- df_cls_centers$cluster_label
  }else{
    print('Using default clone labels in visualization')
    df_cls_centers$cls_lb <- df_cls_centers$cluster_label
  }
  if(col_plt=='cluster_label' | grepl('Lineage',col_plt)){
    cls_col <- 'black'
  }else{
    cls_col <- 'red'
  }
  patient_names <- data.frame(datatag=c('SA609','SA535','SA1035'),
                              pt_names=c('Pt4','Pt5','Pt6'))
  rownames(patient_names) <- patient_names$datatag
  my_font <- "Helvetica"
  plttitle <- patient_names[datatag,'pt_names']
  
  
  p <- ggplot(rd, aes_string(x=x_plt, y=y_plt, color=col_plt)) +
    geom_point(size=1.1, alpha=0.9, shape=21) + #size=0.4
    # viridis::scale_color_viridis(discrete=F, alpha=0.8) +
    scale_color_manual(values = cols_use) +
    theme(plot.title = element_text(color="black", size=17, family=my_font),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text = element_blank(),
          # axis.title = element_text(color="black", size=11),
          axis.title = element_blank(),
          legend.position = 'none',
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white",colour = "white")) #+ 
  # labs(title = plttitle)
  p <- p + geom_curve(
    aes(x = x1, y = y1, xend = x2, yend = y2),
    data = curves, inherit.aes = F,
    curvature = 0.25, size = 0.6, colour = "#2F4F4F",alpha=0.8,
    arrow = arrow(length = unit(0.04, "npc"))
  )
  # Draw cluster labels
  # p <- p + annotate("text", x = df_cls_centers$x1, y = df_cls_centers$y1, 
  #                   label = df_cls_centers$cls_lb, size=6, color=cls_col)
  
  # Draw cluster centers
  for(startp in ls_startp){
    p <- p + annotate(geom="point", startp$x1, y = startp$y1, 
                      colour = "green", size = 6.5, alpha=0.3) # starting point
  }
  for(endp in ls_endp){
    p <- p + annotate(geom="point", endp$x1, y = endp$y1, 
                      colour = "red", size = 6, alpha=0.3) # end point
  }
  
  p <- p + annotate(geom="point", df_cls_centers1$x1, y = df_cls_centers1$y1, 
                    colour = "#2F4F4F", size = 2.5, alpha=0.6) 
  # 
  # p
  
  # png(paste0(output_dir,"ts_slingshot_out_wholedataset_",datatag,".png"), height = 2*400, width=2*600,res = 2*72)
  # print(p)
  # dev.off()
  # ggsave(paste0(save_figs_dir,"umaps_lineages_",datatag,".png"),
  #        plot = p,
  #        height = 4,
  #        width = 5.5,
  #        # useDingbats=F,
  #        dpi=250)
  
  saveRDS(p, paste0(output_dir,"ts_slingshot_out_wholedataset_clones_",datatag,".rds"))
  return(p)
}


plot_whole_lineages_treatmentSt <- function(pt, rd, 
                                           curves, df_cls_centers, 
                                           output_dir, datatag, 
                                           downsample_ratio=0.25, cols_use=NULL,
                                           x_plt='UMAP_1', y_plt='UMAP_2'){
  col_plt <- 'treatmentSt'
  # pt <- pseudo_out 
  # rd_backup <- rd 
  # curves, df_cls_centers, 
  # save_figs_dir, datatag, cols_use=NULL,
  # x_plt='UMAP_1'
  # y_plt='UMAP_2'
  # head(curves)
  if(is.null(cols_use)){
    # res2 <- set_color_treatmentSt(rd[,'treatmentSt'], datatag, cols=NULL)
    # cols_use <- res2$clone_palette
    
    ## Version with time point colors
    # ts_colors <- set_color_treatmentSt(rd[,'treatmentSt'], datatag, cols=NULL)
    # cols_use <- ts_colors$color_code
    # names(cols_use) <- ts_colors$treatmentSt
    
    ## Very simple version
    cols_use <- c("#4d4d4d","#3575e8","#ffd700")
    names(cols_use) <- c('UnRx','Rx','RxH')
  }
  cols_use <- c(cols_use, 'NA'='#E8E8E8')
  
  ls_startp <- list()
  ls_endp <- list()
  for(lid in unique(curves$ligneage_id)){
    curves_tmp <- curves %>% 
      dplyr::filter(ligneage_id==lid)
    # print(dim(curves_tmp))
    obs_cls <- unique(c(curves_tmp$cls1, curves_tmp$cls2))
    df_cls_centers_tmp <- df_cls_centers %>% 
      dplyr::filter(cluster_label %in% obs_cls)
    # class(pt) dim(pt) colnames(pt)
    # print(dim(df_cls_centers_tmp))
    startp <- df_cls_centers_tmp %>% 
      dplyr::filter(cluster_label == curves_tmp$cls1[1])
    ls_startp[[curves_tmp$cls1[1]]] <- startp
    endp <- df_cls_centers_tmp %>% 
      dplyr::filter(cluster_label == curves_tmp$cls2[dim(curves_tmp)[1]])
    ls_endp[[curves_tmp$cls2[dim(curves_tmp)[1]]]] <- endp
    if(!lid %in% colnames(pt)){
      stop('Do not exist any results corresponding to given lineage, double check input data')
    }
  }
  
  # head(df_cls_centers)
  
  
  # View(head(pt))
  # head(rd)
  
  pt$cell_id <- rownames(pt)
  if(!'cell_id' %in% colnames(rd)){
    rd$cell_id <- rownames(rd)
  }
  if(downsample_ratio>0){
    print('Before downsampling')
    print(dim(rd))
    ## Original version
    # cells_use <- rd %>%
    #   # dplyr::filter(cluster_label %in% obs_cls)%>%
    #   dplyr::pull(cell_id)
    
    ## Downsample data, make plot clearer
    cells_use <- downsample_by_cluster_label_treatmentSt(rd, downsample_ratio=downsample_ratio) #0.15
    print(length(cells_use))
    pt <- pt %>%
      dplyr::filter(cell_id %in% cells_use)
  }
  
  # rd <- rd %>% left_join(pt, by=c('cell_id'))
  rd <- rd %>% inner_join(pt, by=c('cell_id'))
  # print('After downsampling')
  print(dim(rd))
  # colors <- pal[cut(rd[,lid], breaks = 100)]
  # rd$treatmentSt <- ifelse(!is.na(rd[,lid]),rd$treatmentSt,'NA')
  
  # length(colors)
  df_cls_centers1 <- df_cls_centers # clusters center
  # Just to avoid overlap with cluster centers
  # df_cls_centers$x1 <- df_cls_centers$x1 + runif(dim(df_cls_centers)[1], -0.5, 0.5)
  # df_cls_centers$y1 <- df_cls_centers$y1 + runif(dim(df_cls_centers)[1], -0.5, 0.5)
  
  # Just to avoid overlap with cluster centers
  df_cls_centers$x1 <- df_cls_centers$x1 + 0.5
  df_cls_centers$y1 <- df_cls_centers$y1 - 0.2
  
  # Just for easy observe output in visualization
  if(datatag=='SA609'){
    # unique(df_cls_centers$cluster_label)
    meta_cluster <- data.frame(cluster_label=c('10','6','1','7','9','4','8','0','5','3','2'),
                                      cls_lb=c('0','1','2','3','4','5','6','7','8','9','10'))
    df_cls_centers <- df_cls_centers %>% inner_join(meta_cluster, by=c('cluster_label'))
  }
  else if(datatag=='SA535'){
    meta_cluster <- data.frame(cluster_label=c('0','1','2','3','4','5','6','7','8','9','10'),
                               cls_lb=c('8','5','2','7','9','6','10','0','4','3','1'))
    df_cls_centers <- df_cls_centers %>% inner_join(meta_cluster, by=c('cluster_label'))
  }else if(datatag=='SA1035'){
    # meta_cluster <- data.frame(cluster_label=c('0','1','2','3','4','5','6','7','8'),
    #                            cls_lb=c('1','0','6','5','3','8','4','7','2')) # version 1 
    meta_cluster <- data.frame(cluster_label=c('0','1','10','6','3','7','9','5','11','4','2','8'),
                               cls_lb=c('0','1','2','3','4','5','6','7','8','9','10','11'))
    df_cls_centers <- df_cls_centers %>% inner_join(meta_cluster, by=c('cluster_label'))
    # df_cls_centers$cls_lb <- df_cls_centers$cluster_label
  }else{
    print('Using default clone labels in visualization')
    df_cls_centers$cls_lb <- df_cls_centers$cluster_label
  }
  if(col_plt=='cluster_label' | grepl('Lineage',col_plt)){
    cls_col <- 'black'
  }else{
    cls_col <- 'red'
  }
  patient_names <- data.frame(datatag=c('SA609','SA535','SA1035'),
                              pt_names=c('Pt4','Pt5','Pt6'))
  rownames(patient_names) <- patient_names$datatag
  my_font <- "Helvetica"
  plttitle <- patient_names[datatag,'pt_names']
  
  if(col_plt=='treatmentSt'){ # simple treatment status version
    rd <- rd %>%
      dplyr::mutate(
        treatmentSt=case_when(
          grepl('T$',treatmentSt) ~ 'Rx',
          grepl('TU$',treatmentSt) ~ 'RxH',
                             TRUE ~ 'UnRx'
        ))
  }
  p <- ggplot(rd, aes_string(x=x_plt, y=y_plt, color=col_plt)) +
    geom_point(size=1.1, alpha=0.9, shape=21) + #size=0.4
    # viridis::scale_color_viridis(discrete=F, alpha=0.8) +
    scale_color_manual(values = cols_use) +
    theme(plot.title = element_text(color="black", size=17, family=my_font),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text = element_blank(),
          # axis.title = element_text(color="black", size=11),
          axis.title = element_blank(),
          legend.position = 'none',
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white",colour = "white")) #+ 
          # labs(title = plttitle)
  p <- p + geom_curve(
    aes(x = x1, y = y1, xend = x2, yend = y2),
    data = curves, inherit.aes = F,
    curvature = 0.25, size = 0.6, colour = "#2F4F4F",alpha=0.8,
    arrow = arrow(length = unit(0.04, "npc"))
  )
  # Draw cluster labels
  # p <- p + annotate("text", x = df_cls_centers$x1, y = df_cls_centers$y1, 
  #                   label = df_cls_centers$cls_lb, size=6, color=cls_col)
  
  # Draw cluster centers
  for(startp in ls_startp){
    p <- p + annotate(geom="point", startp$x1, y = startp$y1, 
                      colour = "green", size = 6.5, alpha=0.3) # starting point
  }
  for(endp in ls_endp){
    p <- p + annotate(geom="point", endp$x1, y = endp$y1, 
                      colour = "red", size = 6, alpha=0.3) # end point
  }
  
  p <- p + annotate(geom="point", df_cls_centers1$x1, y = df_cls_centers1$y1, 
                    colour = "#2F4F4F", size = 2.5, alpha=0.6) 
  # 
  # p
  
  # png(paste0(output_dir,"ts_slingshot_out_wholedataset_",datatag,".png"), height = 2*400, width=2*600,res = 2*72)
  # print(p)
  # dev.off()
  # ggsave(paste0(save_figs_dir,"umaps_lineages_",datatag,".png"),
  #        plot = p,
  #        height = 4,
  #        width = 5.5,
  #        # useDingbats=F,
  #        dpi=250)
  
  saveRDS(p, paste0(output_dir,"ts_slingshot_out_wholedataset_",datatag,".rds"))
  return(p)
}

# pt <- pseudo_out

downsample_by_cluster_label_treatmentSt <- function(rd, downsample_ratio=0.3){
  if(!'cell_id' %in% colnames(rd)){
    rd$cell_id <- rownames(rd)
  }
  set.seed(243)
  cls <- unique(rd$cluster_label)
  cells_use <- c()
  for(cl in cls){
    tmp <- rd %>%
      dplyr::filter(cluster_label==cl) 
    dim(tmp)
    for(s in unique(tmp$treatmentSt)){
      cells_cls <- tmp %>%
        dplyr::filter(treatmentSt==s)  %>%
        dplyr::pull(cell_id)
      if(length(cells_cls)>=10){ # very small outliers, remove from downsampling list
        cells_cls <- sample(cells_cls, round(length(cells_cls)*downsample_ratio))
        cells_use <- c(cells_use, cells_cls)  
      }
      
    }  
  }
  print('Number of downsampling cells: ')
  print(length(cells_use))
  # rd <- rd %>%
  #   dplyr::filter(cell_id %in% cells_use) 
  return(cells_use)
}
plot_given_lineage_treatmentSt <- function(pt, lid, rd, 
                                           curves, df_cls_centers, 
                                           output_dir, datatag, plttitle, cols_use=NULL,
                                           x_plt='UMAP_1', y_plt='UMAP_2', downsample_viz=FALSE){
  # head(curves)
  if(is.null(cols_use)){
    # ts_colors <- set_color_treatmentSt(rd[,'treatmentSt'], datatag, cols=NULL)
    ts_colors <- set_color_treatmentSt_simplify_version(rd[,'treatmentSt'], datatag, cols=NULL)
    cols_use <- ts_colors$color_code
    # names(cols_use) <- ts_colors$label
    # lg2 <- plot_colors(cols_use, save_dir, legend_label='treatment', 1)
    names(cols_use) <- ts_colors$treatmentSt
    # res2 <- set_color_treatmentSt(rd[,'treatmentSt'], datatag, cols=NULL)
    # cols_use <- res2$clone_palette
  }
  cols_use <- c(cols_use, 'NA'='#E8E8E8')
  curves <- curves %>% 
    dplyr::filter(ligneage_id==lid)
  # print(dim(curves))
  # head(df_cls_centers)
  obs_cls <- unique(c(curves$cls1, curves$cls2))
  df_cls_centers <- df_cls_centers %>% 
    dplyr::filter(cluster_label %in% obs_cls)
  # class(pt) dim(pt) colnames(pt)
  # print(dim(df_cls_centers))
  startp <- df_cls_centers %>% 
    dplyr::filter(cluster_label == curves$cls1[1])
  endp <- df_cls_centers %>% 
    dplyr::filter(cluster_label == curves$cls2[dim(curves)[1]])
  if(lid %in% colnames(pt)){
    pt <- pt[,lid, drop=F]
  }else{
    stop('Do not exist any results corresponding to given lineage, double check input data')
  }
  # View(head(pt))
  # head(rd)
  
  pt$cell_id <- rownames(pt)
  if(!'cell_id' %in% colnames(rd)){
    rd$cell_id <- rownames(rd)
  }
  rd <- downsample_by_cluster_label(rd, downsample_ratio=0.4)
  cells_use <- rd %>%
    dplyr::filter(cluster_label %in% obs_cls)%>%
    dplyr::pull(cell_id)
  pt <- pt %>%
    dplyr::filter(cell_id %in% cells_use)
  
  rd <- rd %>% left_join(pt, by=c('cell_id'))
  # dim(rd)
  # colors <- pal[cut(rd[,lid], breaks = 100)]
  rd$treatmentSt <- ifelse(!is.na(rd[,lid]),rd$treatmentSt,'NA')
  col_plt <- 'treatmentSt'
  # length(colors)
  df_cls_centers1 <- df_cls_centers # clusters center
  # Just to avoid overlap with cluster centers
  # df_cls_centers$x1 <- df_cls_centers$x1 + runif(dim(df_cls_centers)[1], -0.5, 0.5)
  # df_cls_centers$y1 <- df_cls_centers$y1 + runif(dim(df_cls_centers)[1], -0.5, 0.5)
  # Just to avoid overlap with cluster centers
  df_cls_centers$x1 <- df_cls_centers$x1 + 0.5
  df_cls_centers$y1 <- df_cls_centers$y1 - 0.2
  
  # Just for easy observe output in visualization
  if(datatag=='SA609'){
    # unique(df_cls_centers$cluster_label)
    meta_cluster <- data.frame(cluster_label=c('10','6','1','7','9','4','8','0','5','3','2'),
                               cls_lb=c('0','1','2','3','4','5','6','7','8','9','10'))
    df_cls_centers <- df_cls_centers %>% inner_join(meta_cluster, by=c('cluster_label'))
  } else if(datatag=='SA535'){
    meta_cluster <- data.frame(cluster_label=c('0','1','2','3','4','5','6','7','8','9','10'),
                               cls_lb=c('8','5','2','7','9','6','10','0','4','3','1'))
    df_cls_centers <- df_cls_centers %>% inner_join(meta_cluster, by=c('cluster_label'))
  }else if(datatag=='SA1035'){
    # meta_cluster <- data.frame(cluster_label=c('0','1','2','3','4','5','6','7','8'),
    #                            cls_lb=c('1','0','6','5','3','8','4','7','2'))  # version 1
    # df_cls_centers <- df_cls_centers %>% inner_join(meta_cluster, by=c('cluster_label'))
    meta_cluster <- data.frame(cluster_label=c('0','1','10','6','3','7','9','5','11','4','2','8'),
                               cls_lb=c('0','1','2','3','4','5','6','7','8','9','10','11'))
    df_cls_centers <- df_cls_centers %>% inner_join(meta_cluster, by=c('cluster_label'))
    # df_cls_centers$cls_lb <- df_cls_centers$cluster_label
  }else{
    print('Using default clone labels in visualization')
    df_cls_centers$cls_lb <- df_cls_centers$cluster_label
  }
  
  if(col_plt=='cluster_label' | grepl('Lineage',col_plt)){
    cls_col <- 'black'
  }else{
    cls_col <- 'red'
  }
  my_font <- "Helvetica"
  plttitle <- NULL
  unique(rd$treatmentSt)
  p <- ggplot(rd, aes_string(x=x_plt, y=y_plt, color=col_plt)) +
    geom_point(size=0.3, alpha=1, shape=21) +
    # viridis::scale_color_viridis(discrete=F, alpha=0.8) +
    scale_color_manual(values = cols_use) +
    theme(#plot.title = element_blank(),
          plot.title = element_text(color="black", size=14, hjust = 0.5, family=my_font),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text = element_blank(),
          # axis.title = element_text(color="black", size=10, family=my_font),
          axis.title = element_blank(),
          legend.position = 'none',
          # axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white",colour = "grey"),
          panel.spacing = unit(c(0, 0, 0, 0), "null"),
          plot.margin = unit(c(0, 0, 0, 0), "null"),
          legend.margin=margin(0,0,0,0),
          # legend.box.margin=margin(-2,-2,-2,-2),
          panel.border = element_rect(colour = "darkgray", fill=NA, size=0.1)
    ) #+ 
    # labs(title = plttitle, x=NULL, y=NULL)
  p <- p + geom_curve(
    aes(x = x1, y = y1, xend = x2, yend = y2),
    data = curves, inherit.aes = F,
    curvature = 0.25, size = 0.7,colour = "#088F8F", alpha=1,
    arrow = arrow(length = unit(0.06, "npc"))
  )
  
  # Draw cluster labels
  # p <- p + annotate("text", x = df_cls_centers$x1, y = df_cls_centers$y1, 
  #                   label = df_cls_centers$cls_lb, size=5, color=cls_col)
  
  # Draw cluster centers
  p <- p + annotate(geom="point", startp$x1, y = startp$y1, 
                    colour = "green", size = 6, alpha=0.4) # starting point
  p <- p + annotate(geom="point", endp$x1, y = endp$y1, 
                    colour = "red", size = 6, alpha=0.4) # end point
  p <- p + annotate(geom="point", df_cls_centers1$x1, y = df_cls_centers1$y1, 
                    colour = "#2F4F4F", size = 4.5, alpha=0.4) 
  # 
  # p
  
  png(paste0(output_dir,"ts_slingshot_out_",lid,"_",datatag,".png"), height = 2*400, width=2*600,res = 2*72)
  print(p)
  dev.off()
  saveRDS(p, paste0(output_dir,"ts_slingshot_out_",lid,"_",datatag,".rds"))
  return(p)
}

# rd: data frame reduction dimensions
# curves: data frame
# rd <- umaps
# curves <- lg_df
# sds <- crv_umap_embed
# lid <- 'Lineage7'
# pt <- slingPseudotime(sds) %>% as.data.frame()
# pt <- pseudo_out

plot_given_lineage <- function(pt, lid, rd, curves, df_cls_centers, 
                               output_dir, x_plt='UMAP_1', y_plt='UMAP_2'){
  # head(curves)
  
  curves <- curves %>% 
    dplyr::filter(ligneage_id==lid)
  # print(dim(curves))
  # View(curves)
  # head(df_cls_centers)
  obs_cls <- unique(c(curves$cls1, curves$cls2))
  df_cls_centers <- df_cls_centers %>% 
    dplyr::filter(cluster_label %in% obs_cls)
  # class(pt) dim(pt) colnames(pt)
  # print(dim(df_cls_centers))
  startp <- df_cls_centers %>% 
    dplyr::filter(cluster_label == curves$cls1[1])
  endp <- df_cls_centers %>% 
    dplyr::filter(cluster_label == curves$cls2[dim(curves)[1]])
  
  if(lid %in% colnames(pt)){
    pt <- pt[,lid, drop=F]
  }else{
    stop('Do not exist any results corresponding to given lineage, double check input data')
  }
  # View(head(pt))
  # head(rd)
  pt$cell_id <- rownames(pt)
  if(!'cell_id' %in% colnames(rd)){
    rd$cell_id <- rownames(rd)
  }
  cells_use <- rd %>%
    dplyr::filter(cluster_label %in% obs_cls)%>%
    dplyr::pull(cell_id)
  pt <- pt %>%
    dplyr::filter(cell_id %in% cells_use)
  rd <- rd %>% left_join(pt, by=c('cell_id'))
  # dim(rd)
  # colors <- pal[cut(rd[,lineage_id], breaks = 100)]
  
  col_plt <- lid
  # length(colors)
  df_cls_centers1 <- df_cls_centers # clusters center
  # Just to avoid overlap with cluster centers
  # df_cls_centers$x1 <- df_cls_centers$x1 + runif(dim(df_cls_centers)[1], -0.5, 0.5)
  # df_cls_centers$y1 <- df_cls_centers$y1 + runif(dim(df_cls_centers)[1], -0.5, 0.5)
  df_cls_centers$x1 <- df_cls_centers$x1 + 0.5
  df_cls_centers$y1 <- df_cls_centers$y1 - 0.2
  cls_col <- 'red'
  if(datatag=='SA609'){
    meta_cluster <- data.frame(cluster_label=c('10','6','1','7','9','4','8','0','5','3','2'),
                               cls_lb=c('0','1','2','3','4','5','6','7','8','9','10'))
    
  }else if(datatag=='SA1035'){
    # meta_cluster <- data.frame(cluster_label=c('0','1','2','3','4','5','6','7','8'),
    #                            cls_lb=c('1','0','6','5','3','8','4','7','2'))  #v1
    meta_cluster <- data.frame(cluster_label=c('0','1','10','6','3','7','9','5','11','4','2','8'),
                               cls_lb=c('0','1','2','3','4','5','6','7','8','9','10','11'))
  }else if(datatag=='SA535'){
    meta_cluster <- data.frame(cluster_label=c('0','1','2','3','4','5','6','7','8','9','10'),
                               cls_lb=c('8','5','2','7','9','6','10','0','4','3','1'))
  }else{
    stop('Need meta cluster info here!!!')
  }
  # meta_cluster <- data.frame(cluster_label=c('0','1','10','6','3','7','9','5','11','4','2','8'),
  #                            cls_lb=c('0','1','2','3','4','5','6','7','8','9','10','11'))
  df_cls_centers <- df_cls_centers %>% inner_join(meta_cluster, by=c('cluster_label'))
  
  # if(col_plt=='cluster_label' | grepl('Lineage',col_plt)){
  #   cls_col <- 'black'
  # }else{
  #   cls_col <- 'red'
  # }
  
  p <- ggplot(rd, aes_string(x=x_plt, y=y_plt, color=col_plt)) +
    geom_point(size=0.4, alpha=0.6, shape=21) +
    viridis::scale_color_viridis(discrete=F, alpha=0.8) +
    # scale_color_manual(values = cols_use) + 
    theme(plot.title = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text = element_blank(),
          axis.title = element_text(color="black", size=12),
          legend.position = 'none',
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white",colour = "white"))
  p <- p + geom_curve(
    aes(x = x1, y = y1, xend = x2, yend = y2),
    data = curves, inherit.aes = F,
    curvature = 0.25, size = 0.8, colour = "#2F4F4F", alpha=0.8,
    arrow = arrow(length = unit(0.07, "npc")) #length = unit(2, "mm")
  )
  # p <- p + annotate("text", x = df_cls_centers$x1, y = df_cls_centers$y1, 
  #                   label = df_cls_centers$cls_lb, size=7, color=cls_col)
  
  # Draw cluster centers
  p <- p + annotate(geom="point", startp$x1, y = startp$y1, 
                    colour = "green", size = 4, alpha=0.4) # starting point
  p <- p + annotate(geom="point", endp$x1, y = endp$y1, 
                    colour = "red", size = 4, alpha=0.4) # end point
  p <- p + annotate(geom="point", df_cls_centers1$x1, y = df_cls_centers1$y1, 
                    colour = "#2F4F4F", size = 2, alpha=0.7) 
  
  # 
  # p
  
  png(paste0(output_dir,"slingshot_out_lineage_",lid,"_",datatag,".png"), height = 2*400, width=2*600,res = 2*72)
  print(p)
  dev.off()
  return(p)
}


# input data: crv1: pseudotime ordering object
# crv <- crv_umap_embed
get_clusters_centers <- function(crv){
  x <- SlingshotDataSet(crv)
  X <- reducedDim(x)
  clusterLabels <- slingClusterLabels(x)
  connectivity <- slingMST(x)
  clusters <- rownames(connectivity)
  nclus <- nrow(connectivity)
  centers <- t(vapply(clusters,function(clID){
    w <- clusterLabels[,clID]
    return(apply(X, 2, weighted.mean, w = w))
  }, rep(0,ncol(X))))
  rownames(centers) <- clusters
  centers_df <- data.frame(centers)
  centers_df$cluster_label <- rownames(centers)
  return(centers_df)
}
# crv1: pseudotime ordering object
# umaps <- reducedDims(sce)[['UMAP']]
# return a dataframe, lineages, and curves: edges x1, y1, x2, y2


get_lineages <- function(crv1, crv_umap_embed=NULL, umaps=NULL){
  if(is.null(crv_umap_embed) & !is.null(umaps)){
    crv_umap_embed <- slingshot::embedCurves(crv1, umaps)
  }
  lgs <- slingLineages(crv1)  
  lg_df <- tibble::tibble()
  for(lg in names(lgs)){
    vtxs <-  lgs[[lg]] 
    for(i in rep(1:(length(vtxs)-1),1)){
      lg_df <- dplyr::bind_rows(lg_df, c(cls1=vtxs[i], cls2=vtxs[i+1], ligneage_id=lg))
    }
  }
  lg_df <- lg_df %>% as.data.frame()
  curves <- get_clusters_centers(crv_umap_embed)
  print(head(curves))
  print(dim(curves))
  
  # Edges: connect one cluster to another cluster
  # Adding center coordinates for first cluster and second cluster
  curves1 <- curves
  colnames(curves1) <- c('x1','y1','cluster_label') 
  curves2 <- curves
  colnames(curves2) <- c('x2','y2','cluster_label') 
  
  lg_df <- lg_df %>% inner_join(curves1, by=c('cls1'='cluster_label'))  
  lg_df <- lg_df %>% inner_join(curves2, by=c('cls2'='cluster_label'))  
  
  print(dim(lg_df))
  print(head(lg_df))
  # df_annot <- curves1
  # df_annot1 <- curves1
  # lg_df <- ls_lineages$lineages
  lg_df$idx <- paste0(lg_df$cls1,'_',lg_df$cls2) # just to plot all curves for whole dataset
  lg_df1 <- lg_df[!duplicated(lg_df$idx),]
  print(dim(lg_df1))
  return(list(lineages=lg_df, curves=curves1, unique_curves=lg_df1))
}



 
# Remove root cluster that is shared by several lineages 
# Get list of cells in each lineage
# Get list of clones, remove clone that less than 2.5% of population 
# Results in lineage and clone list
# If variation > 0 in cis genes, else trans genes 
# Each gene, cis trans lineage 
# Add 3 annotation bars 

get_cis_trans_genes <- function(sce, crv_umap_embed, 
                                output_dir, datatag, 
                                segment_CNV_fn, sign_genes_fn,
                                start_cls=NULL){
  output_dir <- paste0(output_dir,'cis_trans_lineages/')
  if(!dir.exists(output_dir)){
    dir.create(output_dir)  
  }
  
  umaps <- reducedDims(sce)[['UMAP']] %>% as.data.frame()
  umaps$cluster_label <- sce$cluster_label
  umaps$cluster_label <- as.factor(umaps$cluster_label)
  umaps$clone <- get_unique_clone_id(sce$clone)
  
  # TO DO: get clone label in case of double clones labels 
  
  # print(summary(as.factor(umaps$cluster_label)))
  print(summary(as.factor(umaps$clone)))
  none_clones <- c('None','unassigned','Un','Unassigned')
  umaps <- umaps %>%
    tibble::rownames_to_column('cell_id') %>%
    dplyr::filter(!clone %in% none_clones) 
  # head(umaps)
  # Remove root cluster that is shared by several lineages 
  # start_cls <- '10'
  if(!is.null(start_cls)){
    cells_use <- umaps %>%
      # tibble::rownames_to_column('cell_id') %>%
      dplyr::filter(cluster_label!=start_cls) %>%
      dplyr::pull(cell_id)
  }else{
    cells_use <- umaps %>%
      dplyr::pull(cell_id)
  }
  print(length(cells_use))
  # dim(umaps)
  umaps$treatmentSt <- colData(sce)[rownames(umaps),'treatmentSt']
  
  pseudo_out <- slingshot::slingPseudotime(crv_umap_embed) %>% as.data.frame()
  dim(pseudo_out)  
  pseudo_out <- pseudo_out[cells_use,]
  lg_ls <- colnames(pseudo_out)
  
  meta_clone_lg <- tibble::tibble()
  thrs_small_clone <- 0.025 
  for(lg in colnames(pseudo_out)){
    cells_lg <- pseudo_out %>%
      tibble::rownames_to_column('cell_id') %>%
      dplyr::select(!!sym(lg), cell_id) %>%
      dplyr::filter(!is.na(!!sym(lg))) %>%
      dplyr::pull(cell_id)
    nb_cells <- length(cells_lg)
    stat_clones_lg <- umaps %>%
      dplyr::filter(cell_id %in% cells_lg) %>%
      dplyr::group_by(clone) %>%
      dplyr::summarise(pct_cells=n()/nb_cells) %>%
      dplyr::filter(pct_cells>thrs_small_clone) %>%
      dplyr::pull(clone)
    tmp <- tibble::tibble(lineage=rep(lg, length(stat_clones_lg)),clone=stat_clones_lg)
    meta_clone_lg <- dplyr::bind_rows(meta_clone_lg, tmp)
  }
  
  # Get lineage description for plotting
  plttlts <- get_lineage_plot_label(datatag)
  meta_lineage <- tibble::tibble(lineage=paste0('Lineage',rep(1:length(plttlts),1)),
                                 lineage_desc=plttlts)
  meta_clone_lg <- meta_clone_lg %>% left_join(meta_lineage, by='lineage')
  data.table::fwrite(meta_clone_lg, paste0(output_dir,datatag,'_meta_clone_lineages.csv'))
  # meta_clone_lg <- data.table::fread(paste0(output_dir,datatag,'_meta_clone_lineages.csv')) %>% as.data.frame()
  # Get segment CNV profile
  cnv <- data.table::fread(segment_CNV_fn) %>% as.data.frame()
  cnv$ens_gene_id <- cnv$V1
  cnv$V1 <- NULL
  # dim(cnv)
  # head(cnv)
  
  signf_genes <- data.table::fread(sign_genes_fn) %>% as.data.frame()
  common_genes <- intersect(signf_genes$ens_gene_id, cnv$ensembl_gene_id)
  
  print('Number of common genes: ')
  print(length(common_genes))
  dim(signf_genes)
  genes_stat <- tibble::tibble()
  # 'Unmapped CNV'
  for(lg in lg_ls){
    obs_clones <- meta_clone_lg %>%
      dplyr::filter(lineage==lg) %>%
      dplyr::pull(clone)
    
    for(g in common_genes){
      cnv_profile <- cnv %>%
        dplyr::filter(ensembl_gene_id==g) %>%
        dplyr::select(all_of(obs_clones))
      cnv_profile <- unlist(cnv_profile[1,])
      cnv_profile <- cnv_profile[!is.na(cnv_profile)]
      if(!is.null(cnv_profile) & length(cnv_profile)>1){
        gene_type <- ifelse(var(cnv_profile)>0,'cis','trans')  #c: cis, t: trans
      }else{
        gene_type <- 'trans'
      }
      
      tmp <- c(ens_gene_id=g, lineage=lg, gene_type=gene_type)
      genes_stat <- dplyr::bind_rows(genes_stat, tmp)
    }
  }
  # genes_stat$gene_type <- ifelse(is.na(genes_stat$gene_type),'trans',genes_stat$gene_type)
  dim(genes_stat)
  unmapped_genes <- signf_genes$ens_gene_id[!signf_genes$ens_gene_id %in% common_genes]
  length(unmapped_genes)
  nb_lg <- length(lg_ls)
  
  for(g in unmapped_genes){
    tmp <- tibble::tibble(ens_gene_id=rep(g,nb_lg), lineage=lg_ls, gene_type=rep('unmapped',nb_lg))
    genes_stat <- dplyr::bind_rows(genes_stat, tmp)
  }  
    
  print(summary(as.factor(genes_stat$gene_type)))
  data.table::fwrite(genes_stat, paste0(output_dir,datatag,'_genes_cis_trans_lineages.csv'))
  
}


plot_all_lingeages <- function(sce, crv_umap_embed, output_dir, datatag){
  # output_dir <- "/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/slingshot_output_400/"
  # dir.create(output_dir)
  # save_figs_dir <- paste0(output_dir,'figs_v4/')
  # save_figs_dir <- paste0(output_dir,'figs_v5/')
  save_figs_dir <- paste0(output_dir,'figs_revision/')
  if(!dir.exists(save_figs_dir)){
    dir.create(save_figs_dir)
  }  
  ls_lineages <- get_lineages(crv_umap_embed, crv_umap_embed, umaps=NULL)
  unique_curves_ls <- ls_lineages$unique_curves
  lg_df <- ls_lineages$lineages
  
  umaps <- reducedDims(sce)[['UMAP']] %>% as.data.frame()
  umaps$cluster_label <- sce$cluster_label
  umaps$cluster_label <- as.factor(umaps$cluster_label)
  # umaps$treatmentSt <- colData(sce)[rownames(umaps),'treatmentSt']
  umaps$treatmentSt <- sce$treatmentSt
  
  # umaps$clone <- sce$clone
  dim(umaps)
  umaps$cell_id <- rownames(umaps)
  umaps$cell_id[1]
  
  
  ## Revision, need unique clonal labels here. 
  ## To Do: get unique labels from here, and change labels in sce file? 
  ## SA609
  clonal_df <- data.table::fread('/home/htran/storage/datasets/drug_resistance/rna_results/SA609_rna/slingshot_trajectory/withBE_SA609_v2/clone_labels_unique_SA609.csv.gz')
  ## SA535
  # clonal_df <- data.table::fread('/home/htran/storage/datasets/drug_resistance/rna_results/SA535_rna/slingshot_trajectory/clone_labels_unique_SA535.csv.gz')

  ## SA1035
  # clonal_df <- data.table::fread('/home/htran/storage/datasets/drug_resistance/rna_results/SA1035_rna/slingshot_trajectory/clone_labels_unique_SA1035.csv.gz')
  # dim(clonal_df)
  
  # colnames(clonal_df)
  # clonal_df$cell_id[1]
  # unique(clonal_df$clone)
  clonal_df <- clonal_df %>%
    dplyr::select(cell_id, clone) %>%
    dplyr::mutate(
      clone = case_when(
        clone=='unassigned' ~ 'None', 
        TRUE ~ clone
      )
    )
  umaps <- umaps %>%
    dplyr::inner_join(clonal_df, by='cell_id')
  dim(umaps)
  unique(umaps$clone)
  
  # Draw clusters
  res <- set_color_clusters(umaps$cluster_label)
  lg <- plot_colors(res$clone_palette, save_figs_dir, legend_label='clusters')
  cols_use <- res$clone_palette
  pc <- plot_curves(datatag, umaps, unique_curves_ls, 
                    ls_lineages$curves, cols_use, save_figs_dir, col_plt='cluster_label', 
                    x_plt='UMAP_1', y_plt='UMAP_2')
  
  # Draw clone label
  # res1 <- set_color_clones(umaps$clone)
  cols_use_clones <- get_color_clones_v2(datatag, unique(umaps$clone))
  # cols_use <- res1$clone_palette
  lg1 <- plot_colors(cols_use_clones, save_figs_dir, legend_label='clone')
  
  pcl <- plot_curves(datatag, umaps, unique_curves_ls, 
                     ls_lineages$curves, cols_use_clones, save_figs_dir, col_plt='clone', 
                     x_plt='UMAP_1', y_plt='UMAP_2')
  
  # Draw treatment
  # sce$treatmentSt <- sce$treat
  ts_colors <- set_color_treatmentSt(colData(sce)[,'treatmentSt'], datatag, cols=NULL)
  ts_colors <- set_color_treatmentSt_simplify_version(colData(sce)[,'treatmentSt'], datatag, cols=NULL)
  # data.table::fwrite(ts_colors, paste0(save_figs_dir, 'treatment_meta_colors.csv'))
  cols_use <- ts_colors$color_code
  names(cols_use) <- ts_colors$label
  lg2 <- plot_colors(cols_use, save_figs_dir, legend_label='treatment', 3)
  
  # ## For simplify version SA609
  # cols_use <- c("#989898","#4d4d4d","#00FFFF","#7aa4f0","#D5FF00","#FFFF00")
  # # names(cols_use) <- c('1 UnRx', '2,3,4 UnRx','1 Rx', '2,3,4 Rx','1 RxH', '2,3,4 RxH')
  # names(cols_use) <- c('UnRx-X4', 'UnRx-X5,X6,X7','Rx-X4', 'Rx-X5,X6,X7','RxH-X5', 'RxH-X6,X7')
  
  ## For simplify-est version SA609 with only 3 status
  cols_use <- c("#4d4d4d","#3575e8","#ffd700")
  names(cols_use) <- c('UnRx','Rx','RxH')
  lg2 <- plot_colors(cols_use, save_figs_dir, legend_label='Treatment', 1)
  
  ## For simplify version SA535
  cols_use <- c("#989898","#4d4d4d","#00FFFF","#7aa4f0","#D5FF00","#FFFF00")
  names(cols_use) <- c('UnRx-X6', 'UnRx-X7,X8,X9','Rx-X6', 'Rx-X7,X8,X10','RxH-X7', 'RxH-X8,X10')

  ## For simplify version SA1035
  # cols_use <- c("#989898","#4d4d4d","#00FFFF","#7aa4f0","#D5FF00","#FFFF00")
  # # names(cols_use) <- c('1 UnRx', '2,3,4 UnRx','1 Rx', '2,3,4 Rx','1 RxH', '2,3,4 RxH')
  # names(cols_use) <- c('UnRx-X4', 'UnRx-X5,X6,X7,X8','Rx-X5', 'Rx-X6,X7,X8','RxH-X6', 'RxH-X7,X8')
  
  ## For simplify-est version SA1035 with only 3 status
  cols_use <- c("#4d4d4d","#3575e8","#ffd700")
  names(cols_use) <- c('UnRx','Rx','RxH')
  lg2 <- plot_colors(cols_use, save_figs_dir, legend_label='Treatment', 1)
  
  # cols_use <- c("#989898","#4d4d4d","#00FFFF","#7aa4f0","#D5FF00","#FFFF00")
  # names(cols_use) <- c('1 UnRx', '2,3,4 UnRx','1 Rx', '2,3,5 Rx','1 RxH', '2,4 RxH')
  lg2 <- plot_colors(cols_use, save_figs_dir, legend_label='', nrow_plt=2) #Treatment\n Cycle
  lg2 <- plot_colors(cols_use, save_figs_dir, legend_label='', 6) #Treatment\n Cycle
  # lg2 <- plot_colors(cols_use, save_figs_dir, legend_label='treatment', round(dim(ts_colors)[1]/4,0))
  
  names(cols_use) <- ts_colors$treatmentSt
  pt <- plot_curves(datatag, umaps, unique_curves_ls, 
                    ls_lineages$curves, cols_use, save_figs_dir, col_plt='treatmentSt', 
                    x_plt='UMAP_1', y_plt='UMAP_2')
  
  
  
  # Draw different lineages
  df_cls_centers <- ls_lineages$curves
  rd <- umaps
  # p <- plot_lingeages(umaps, lg_df1)
  curves <- ls_lineages$lineages
  pseudo_out <- slingPseudotime(crv_umap_embed) %>% as.data.frame()
  # lid <- 'Lineage2'
  
  # plt_ls <- list()
  # for(lid in unique(curves$ligneage_id)){
  #   plt_ls[[lid]] <- plot_given_lineage(pseudo_out, lid, rd, curves, df_cls_centers, 
  #                                       output_dir, x_plt='UMAP_1', y_plt='UMAP_2')
  #   
  # }
  # plt_ls$Lineage1
  
  plttlts <- get_lineage_plot_label(datatag)
  lgs <- unique(curves$ligneage_id)
  
  ## Plotting each lineage
  plt_ts_ls <- list()
  for(i in rep(1:length(lgs),1)){
    lid <- lgs[i]
    if(!is.null(plttlts)){
      plttitle <- plttlts[i]
    }else{
      plttitle <- gsub('Lineage','L',lid)  #Using default lineages name
    }
    
    plt_ts_ls[[lid]] <- plot_given_lineage_treatmentSt(pseudo_out, lid, rd, 
                                                       curves, df_cls_centers, 
                                                       save_figs_dir, datatag, plttitle, 
                                                       cols_use=NULL,
                                                       x_plt='UMAP_1', y_plt='UMAP_2')
    
  }
  
  # plt_ts_ls$Lineage1
  # plg <- cowplot::plot_grid(lg2, NULL, ncol=2, rel_widths = c(1,1.5))
    ## SA535
  # p <- cowplot::plot_grid(plt_ts_ls$Lineage4, plt_ts_ls$Lineage2, 
  #                         plt_ts_ls$Lineage1, plt_ts_ls$Lineage3, nrow=2) + 
  #   theme(plot.background = element_rect(fill = "white", colour = "white"))
  # p
  ## SA1035
  p <- cowplot::plot_grid(plt_ts_ls$Lineage1, plt_ts_ls$Lineage2,
                          plt_ts_ls$Lineage3, plt_ts_ls$Lineage4, nrow=2) +
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  # p
  dim(pseudo_out)
  dim(rd)
  
  p <- plot_whole_lineages_treatmentSt(pseudo_out, rd, 
                                  curves, df_cls_centers, 
                                  save_figs_dir, datatag, 
                                  downsample_ratio=0.3, cols_use=NULL,
                                  x_plt='UMAP_1', y_plt='UMAP_2')
  p <- plot_whole_lineages_clones(pseudo_out, rd, 
                                   curves, df_cls_centers, 
                                   save_figs_dir, datatag, 
                                   downsample_ratio=0.2, cols_use=cols_use_clones,
                                   x_plt='UMAP_1', y_plt='UMAP_2')
  # p
  # p_total <- p + lg1
  # p_total
  ggsave(paste0(save_figs_dir,"SUPP_Fig11_trajectory_SA1035_plg_lineage.svg"),
         plot = p,
         height = 2.5,
         width = 3.5,
         # useDingbats=F,
         dpi=150)
  
  p1 <- cowplot::plot_grid(plt_ts_ls$Lineage3,plt_ts_ls$Lineage2, plt_ts_ls$Lineage1, ncol=1)
  ptotal <- cowplot::plot_grid(p,p1, rel_widths = c(2,1))
  
  # SA535
  # p1 <- cowplot::plot_grid(plt_ts_ls$Lineage4,plt_ts_ls$Lineage2, ncol=1) 
  # p2 <- cowplot::plot_grid(plt_ts_ls$Lineage1, plt_ts_ls$Lineage3, NULL, ncol=3)
  # ptotal1 <- cowplot::plot_grid(p,p1, rel_widths = c(2,1), nrow = 1)
  # ptotal <- cowplot::plot_grid(ptotal1,p2, rel_heights = c(2,1), ncol=1)
  
  png(paste0(output_dir,"slingshot_out_summary_",datatag,".png"), height = 2*500, width=2*700,res = 2*72)
  print(ptotal)
  dev.off()
  
  res <- list(plt_ls=plt_ls, plt_ts_ls=plt_ts_ls, umaps=umaps, ls_lineages=ls_lineages)
  saveRDS(res, paste0(output_dir,'res_slingshot_plots.rds'))
  # return(res)  
}

get_lineage_plot_label <- function(datatag=''){
  plttlts <- NULL
  if(datatag=='SA609'){
    # meta_lineage <- data.frame(lineage=c("Lineage 1","Lineage 2","Lineage 3"),
    #                            lineage_desc=c("Lineage 3: Rx","Lineage 2: RxH","Lineage 1: UnRx"))
    # plttlts <- c('Lineage 3: Rx','Lineage 2: RxH','Lineage 1: UnRx')
    plttlts <- c('L3-Rx','L2-RxH','L1-UnRx')
  }else if(datatag=='SA535'){
    # meta_lineage <- data.frame(lineage=c("Lineage 1","Lineage 2","Lineage 3","Lineage 4"),
    #                            lineage_desc=c("Lineage 3: Rx,RxH","Lineage 2: Rx,RxH",
    #                                           "Lineage 4: Rx,1Rx","Lineage 1: UnRx"))
    plttlts <- c("L3-Rx,RxH","L2-Rx,RxH",
                 "L4-Rx,1Rx","L1-UnRx")
  }else if(datatag=='SA1035'){
    # meta_lineage <- data.frame(lineage=c("Lineage1","Lineage2","Lineage3","Lineage4","Lineage5"),
    #                            lineage_desc=c("L4: Rx,RxH","L3: UnRx(Sen), RxH",
    #                                           "L2: Rx(Res), RxH","L1: 1Rx, Rx", "L5:Rx, UnRx"))
    
    # plttlts <- c("L2: Rx,RxH","L4: RxH, UnRx",
    #              "L1: 1Rx,Rx","L3: Rx")
    plttlts <- c("L4-Rx,RxH","L3-UnRx(Sen),RxH","L2-Rx(Res),RxH","L1-1Rx,Rx","L5-Rx,UnRx")
    # plttlts <- NULL
  }else{
    print('Using default lineages name')
    # plttlts <- NULL
  }  
  # }else if(datatag=='SA535'){
  #   # meta_lineage <- data.frame(lineage=c("Lineage 1","Lineage 2","Lineage 3","Lineage 4"),
  #   #                            lineage_desc=c("Lineage 3: Rx,RxH","Lineage 2: Rx,RxH",
  #   #                                           "Lineage 4: Rx,1Rx","Lineage 1: UnRx"))
  #   plttlts <- c("L3-Rx,RxH","L2-Rx,RxH",
  #                "L4-Rx,1Rx","L1-UnRx")
  # }
  
  return(plttlts)
}

get_cluster_annotation_plt <- function(datatag=''){
  if(datatag=='SA609'){
    meta_cluster <- data.frame(cluster_label=c('10','6','1','7','9','4','8','0','5','3','2'),
                               cls_lb=c('0','1','2','3','4','5','6','7','8','9','10'))
    
  }else if(datatag=='SA1035'){
    # meta_cluster <- data.frame(cluster_label=c('0','1','2','3','4','5','6','7','8'),
    #                            cls_lb=c('1','0','6','5','3','8','4','7','2'))  #v1
    meta_cluster <- data.frame(cluster_label=c('0','1','10','6','3','7','9','5','11','4','2','8'),
                               cls_lb=c('0','1','2','3','4','5','6','7','8','9','10','11'))
  }else if(datatag=='SA535'){
    meta_cluster <- data.frame(cluster_label=c('0','1','2','3','4','5','6','7','8','9','10'),
                               cls_lb=c('8','5','2','7','9','6','10','0','4','3','1'))
  }else{
    stop('Need meta cluster info here!!!')
  }
  return(meta_cluster)
}
# rd: data frame reduction dimensions
# curves: data frame
# rd <- umaps
# curves <- lg_df1
# curves <- unique_curves_ls
# df_annot <- ls_lineages$curves
plot_curves <- function(datatag, rd, curves, df_annot, cols_use, output_dir, col_plt='cluster_label', 
                        x_plt='UMAP_1', y_plt='UMAP_2'){
  df_annot1 <- df_annot # clusters center
  # Just to avoid overlap with cluster centers
  # df_annot$x1 <- df_annot$x1 + runif(dim(df_annot)[1], -0.5, 0.5)
  # df_annot$y1 <- df_annot$y1 + runif(dim(df_annot)[1], -0.5, 0.5)
  df_annot$x1 <- df_annot$x1 + 0.5
  df_annot$y1 <- df_annot$y1 - 0.2
  cls_col <- 'red'
  # if(col_plt=='cluster_label'){
  #   cls_col <- 'black'
  # }else{
  #   cls_col <- 'red'
  # }
  if(dim(rd)[1]>15000){
    cell_size <- 0.15
  }else{
    cell_size <- 0.4
  }
  
  meta_cluster <- get_cluster_annotation_plt(datatag)  
  df_annot <- df_annot %>% inner_join(meta_cluster, by=c('cluster_label'))
  p <- ggplot(rd, aes_string(x=x_plt, y=y_plt, color=col_plt)) +
    geom_point(size=cell_size, alpha=0.6, shape=21) + 
    scale_color_manual(values = cols_use) + 
    theme(plot.title = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text = element_blank(),
          axis.title = element_text(color="black", size=12),
          legend.position = 'bottom',
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white",colour = "white"))
  p <- p + guides(color = guide_legend(nrow = 2, override.aes = list(size=6))) #
  p <- p + geom_curve(
    aes(x = x1, y = y1, xend = x2, yend = y2),
    data = curves, inherit.aes = F,
    curvature = 0.25, size = 0.8,colour = "#2F4F4F", alpha=0.8,
    arrow = arrow(length = unit(0.07, "npc"))
  )
  # Draw curves
  # p <- p + annotate("text", x = df_annot$x1, y = df_annot$y1, 
  #                   label = NULL, size=5, color=cls_col) #label = df_annot$cls_lb
  
  # Draw cluster centers
  p <- p + annotate(geom="point", df_annot1$x1, y = df_annot1$y1, 
                    colour = "black", size = 1.5, alpha=0.5) 
  # p
  png(paste0(output_dir,"slingshot_out_curves_",datatag,"_",col_plt,".png"), height = 2*400, width=2*600,res = 2*72)
  print(p)
  dev.off()
  return(p)
}
get_color_clones_v2 <- function(tag, clones_ls, color_fn=NULL){
  if(is.null(color_fn)){
    base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
    color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v4.csv.gz')  
  }
  color_df <- data.table::fread(color_fn)
  # color_df <- data.table::fread(paste0(output_dir,'colorcode_total.csv'))
  color_df <- color_df %>%
    dplyr::filter(datatag==tag)
  
  if(dim(color_df)[1] > 0){
    cols_use <- color_df$colour
    names(cols_use) <- color_df$clone_id
    if(length(clones_ls)>length(cols_use)){
      print('Double check, lack of defined clones in the color file!!!')
      print('List of clones: ')
      print(clones_ls)
      print('List of defined colors: ')
      print(cols_use)
    }
    clones_ls <- clones_ls[clones_ls %in% names(cols_use)]
    cols_use <- cols_use[clones_ls]
    if(!'None' %in% names(cols_use)){
      cols_use['None'] <- '#D3D3D3'  
    }
    # if('None' %in% names(cols_use)){
    #     
    # }
    
  }else{
    # cols_use <- make_clone_palette(obs_clones)
    stop('Error, check color mapping')
  }
  print("Color setting: ")
  print(cols_use)
  return(cols_use)
}

set_color_clones <- function(clusters_labels){
  clusters_label <- gtools::mixedsort(unique(clusters_labels))
  print(clusters_label)
  # colorcode_fn <- "/home/htran/storage/datasets/metastasis_results/SA919X7_new_encoding/SA919_Tyler_wholedata/config/colorcode.csv"
  base_dir <- '/home/htran/Projects/farhia_project/drug_resistant_material/'
  colorcode_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v4.csv.gz')
  
  
  if(file.exists(colorcode_fn)){
    color_df <- data.table::fread(colorcode_fn)
    color_df <- color_df %>%
      dplyr::filter(nb_clones==length(clusters_label))
    clone_palette <- color_df$color
  }else{
    clone_palette <- RColorBrewer::brewer.pal(9,"Set1")
  }
  
  
  if(length(clusters_label)>length(clone_palette)){
    stop('Number of predefined colors smaller than number of input clusters')
  }
  clone_palette <- clone_palette[1:length(clusters_label)]
  names(clone_palette) <- clusters_label
  cols_use <- clone_palette[clusters_labels]
  return(list(clone_palette=clone_palette,cols_use=cols_use))
  
}

set_color_clusters <- function(clusters_labels){
  clusters_label <- gtools::mixedsort(unique(clusters_labels))
  print(clusters_label)
  clone_palette1 <- c("#be5f72", "#d74058", "#dc4229", "#a6552c", "#df956f", "#e47a33",
                     "#d49f34", "#836e2c", "#b2ad5a", "#92b539", "#4c7d38", "#4dc041",
                     "#5dba7f", "#47b8c3", "#6280ca", "#7b57db", "#ce8bd1", "#934f94",
                     "#cb48cb", "#d74391")
  clone_palette2 <- c("#be5f72", "#d74058", "#6280ca", "#00008b", "#5dba7f", "#e47a33",
                     "#d49f34", "#ce8bd1", "#b2ad5a", "#92b539","#47b8c3", "#7b57db",
                     "#cb48cb","#934f94", "#d74391")
  if(length(clusters_label)<=7){
    clone_palette <- RColorBrewer::brewer.pal(9,"Set1")
  }
  else if(length(clusters_label)<=15 & length(clusters_label)>7){
    clone_palette <- clone_palette2
  }else{
    clone_palette <- clone_palette1
  }
  if(length(clusters_label)>length(clone_palette)){
    stop('Number of predefined colors smaller than number of input clusters')
  }
  clone_palette <- clone_palette[1:length(clusters_label)]
  names(clone_palette) <- clusters_label
  cols_use <- clone_palette[clusters_labels]
  return(list(clone_palette=clone_palette,cols_use=cols_use))
  
}
plot_colors <- function(cols, save_dir, legend_label='', nrow_plt=5){
  # clusters_use <- gtools::mixedsort(names(cols))
  # clusters_use <- c('UU','UUU','UUUU','UUUUU','UUT','UUTT','UUTTT','UUTTTTT','UUTU','UUTTU','UUTTTTU')
  # cols <- cols[clusters_use]
  color_df <- data.frame(color_code=cols, cluster=names(cols), vals=c(rep(1:length(cols),1)))
  data.table::fwrite(color_df, paste0(save_dir, legend_label,'_metacolors.csv'))
  # color_df$cluster <- factor(color_df$cluster, levels = lvs)
  color_df$cluster <- factor(color_df$cluster, levels = color_df$cluster)
  p <- ggplot(color_df, aes(x=cluster, y=vals, color=cluster)) +
    # geom_line(aes_string(color=plottype)) +  #,color=colorcode
    geom_point(size=7.5) +
    scale_color_manual(values = cols) + labs(colour='') + 
    theme(legend.text=element_text(size=11, hjust = 0, family=my_font),
          legend.title=element_text(size=11, hjust = 0.5, family=my_font),
          panel.spacing = unit(c(0, 0, 0, 0), "null"),
          legend.key=element_blank())#panel.margin = unit(c(0, 0, 0, 0), "null")
  p <- p + guides(color = guide_legend(title = legend_label, nrow=nrow_plt, override.aes = list(size=5))) 
  
  lg <- cowplot::get_legend(p)
  pc <- cowplot::ggdraw() + cowplot::draw_plot(lg) #+ theme(plot.background = element_rect(fill = "white", colour = "white"))
  saveRDS(pc, paste0(save_dir, gsub(' ','_',legend_label),'_plt.rds'))
  # saveRDS(lg, paste0(save_dir, legend_label,'_plt_legend.rds'))
  # png(paste(save_dir,"lg_",legend_label,".png",sep=""),height = 2*400, width=2*300,res = 2*72)
  # # print(grid.arrange(grobs = plots[select_grobs(layout)], layout_matrix = layout,
  # #                    bottom=" ",right=" "))
  # print(pc)
  # dev.off()
  # ggsave(paste0(save_dir,"pseudo_gene_modules_",datatag, "_lg_color.png"),
  #        plot = cowplot::ggdraw() + cowplot::draw_plot(lg) + 
  #          theme(plot.background = element_rect(fill = "transparent", colour = "white")),
  #        height = 4.5,
  #        width = 3,
  #        # useDingbats=F,
  #        dpi=150)
  # ggsave(paste0(save_dir,"pseudo_gene_modules_",datatag, "_lg_color.svg"),
  #        plot = cowplot::ggdraw() + cowplot::draw_plot(lg) + 
  #          theme(plot.background = element_rect(fill = "transparent", colour = "white")),
  #        height = 4.5,
  #        width = 3,
  #        # useDingbats=F,
  #        dpi=150)
  # ggsave(paste0(save_figs_dir,"SA609_umaps_legends.svg"),
  #        plot = lg2,
  #        height = 1,
  #        width = 2.5,
  #        # useDingbats=F,
  #        dpi=20)
  
  return(pc)
}
get_start_cluster <- function(meta_info, datatag,root_ts=NULL){
  if(is.null(root_ts)){
    if(datatag=='SA535'){
      root_ts <- 'UU'   
    }else{
      root_ts <- 'U'   
    }  
  }
  print(root_ts)
  start_cls_ls <- c()
  for(r in root_ts){
    plt_feature <- colnames(meta_info)[grepl('treat',colnames(meta_info))]
    # colnames(meta_info)[which(colnames(meta_info)==plt_feature)] <- 'treatmentSt'
    tmp <- table(meta_info[,'cluster_label'], meta_info[,plt_feature]) %>% as.data.frame()
    colnames(tmp) <- c('cluster_label','treatmentSt','Freq')
    tmp <- tmp %>%
      dplyr::filter(treatmentSt==r & Freq > 20)
    print(tmp)
    # tmp <- tmp %>%
    #   dplyr::filter(cluster_label==18 & Freq > 20)
    start_cls <- tmp %>%
      dplyr::filter(Freq == max(Freq)) %>%
      dplyr::pull(cluster_label)
    
    # start_cls <- '12'
    start_cls_ls <- c(start_cls_ls, as.character(start_cls))
  }
  names(start_cls_ls) <- root_ts
  return(start_cls_ls)
}

prepare_data_Seurat <- function(sce, output_dir, datatag, save_srt=FALSE){
  reducedDims(sce) <- NULL
  nb_pcs <- 30
  nfeatures_use <- 3000
  if('ID' %in% colnames(rowData(sce)) & !grepl('ENSG',rownames(sce)[1])){
    print(rowData(sce)$ID[1])
    rownames(sce) <- rowData(sce)$ID
  }
  print(rownames(sce)[1])
  exprs <- "logcounts"
  dim(logcounts(sce))
  dim(counts(sce))
  srt <- Seurat::as.Seurat(sce, counts = "counts", data = "logcounts")
  # srt <- readRDS(paste0(output_dir, datatag,'_',nfeatures_use,'_srt.rds'))
  # srt <- Seurat::FindVariableFeatures(object = srt, verbose = T)
  # srt <- Seurat::ScaleData(object = srt, verbose = FALSE) # problematics
  # srt[["RNA"]]@scale.data <- as.matrix(logcounts(sce))
  srt <- SetAssayData(object = srt, slot = "scale.data", new.data = as.matrix(logcounts(sce)))
  print(dim(srt))
  # srt[["RNA"]]@scale.data <- as.matrix(assay(sce, exprs))
  # class(srt[["RNA"]]@data)
  srt <- Seurat::FindVariableFeatures(object = srt, verbose = FALSE,
                                      selection.method = "vst", nfeatures = nfeatures_use)
  var_genes <- Seurat::VariableFeatures(object = srt)
  print(paste0('Nb var genes: ',length(var_genes)))
  
  
  var_genes_df <- data.frame(var_gene=var_genes)
  write.csv(var_genes_df, file=paste0(output_dir,'var_genes_',nfeatures_use,'.csv'), quote=F, row.names = F)
  
  # sce <- sce[as.character(var_genes),]
  # print(dim(sce))
  # if(!save_srt){
  #   saveRDS(sce, paste0(output_dir, datatag,'_',nfeatures_use,'_sce.rds'))
  # }
  
  # Filter low quality clusters, and re-run PCA, UMAP
  # colnames(meta_info)
  # meta_info <- meta_info %>% 
  #   filter(!seurat_clusters %in% c(5, 7))
  # dim(meta_info)
  # srt <- srt[,meta_info$cell_id]
  # dim(srt)
  print("Run PCA")
  srt <- Seurat::RunPCA(object = srt, verbose = FALSE, npcs = nb_pcs) #, features = var_genes
  dims = 1:nb_pcs
  
  srt <- Seurat::FindNeighbors(srt, dims = dims, verbose = FALSE)
  
  srt <- Seurat::FindClusters(srt, resolution = 0.7) # 4: Leiden algo, here using default: Louvain
  
  meta_info <- srt@meta.data
  meta_info$cell_id <- rownames(meta_info)
  data.table::fwrite(meta_info, paste0(output_dir, datatag,'_',nfeatures_use, "_meta_cells.csv"))
  
  print("Run UMAP")
  # # srt <- RunTSNE(object = srt, verbose = FALSE)
  # umap method uwot
  srt <- Seurat::RunUMAP(object = srt, dims = dims, verbose = FALSE) #umap.method = "umap-learn", metric='correlation'
  dim(srt)
  # if(save_srt){
  #   saveRDS(srt, paste0(output_dir, datatag,'_',nfeatures_use,'_srt.rds'))
  # }
  
  pca_df <- as.data.frame(Seurat::Embeddings(object = srt, reduction = "pca"))
  umap_df <- as.data.frame(Seurat::Embeddings(object = srt, reduction = "umap"))
  # colnames(pca_df)
  
  meta_info1 <- meta_info %>%
    dplyr::select(cell_id, clone, Barcode, library_id, sample, treatmentSt, timepoint) #, clone, treatmentSt, treat
  
  umap_df$cell_id <- rownames(umap_df)
  umap_df <- umap_df %>% inner_join(meta_info1, by=c("cell_id"))
  # umap_df$cell_id <- paste0(umap_df$library_id,'_',umap_df$Barcode)
  
  pca_df$cell_id <- rownames(pca_df)
  pca_df <- pca_df %>% inner_join(meta_info1, by=c("cell_id"))
  # pca_df$cell_id <- paste0(pca_df$library_id,'_',pca_df$Barcode)
  print(dim(umap_df))
  print(dim(pca_df))
  
  data.table::fwrite(pca_df, paste0(output_dir, datatag,'_',nfeatures_use, "_norm_pca.csv"))
  data.table::fwrite(umap_df, paste0(output_dir, datatag,'_',nfeatures_use, "_norm_umap.csv"))
  
  sce <- sce[var_genes,]  
  print(dim(sce))
  # colnames(pca_df) <- ifelse(grepl('pca_', colnames(pca_df)),
  #                            gsub('pca','PC',colnames(pca_df)), colnames(pca_df))
  # 
  # colnames(umap_df) <- ifelse(grepl('umap_', colnames(umap_df)),
  #                            gsub('umap','UMAP',colnames(umap_df)), colnames(umap_df))
  # 
  pcs <- pca_df[,paste0('PC_', rep(1:nb_pcs,1))]
  umaps <- umap_df[,paste0('UMAP_', rep(1:2,1))]
  colnames(srt)[1]
  sce <- sce[,colnames(srt)]
  reducedDims(sce) <- SimpleList(PCA = as.matrix(pcs), UMAP = as.matrix(umaps))
  dim(sce)
  dim(srt)
  # rownames(meta_info) <- meta_info$cell_id
  colData(sce)$cluster_label <- meta_info[colnames(sce),'seurat_clusters']
  print(summary(as.factor(colData(sce)$cluster_label)))
  saveRDS(sce, paste0(output_dir, datatag,'_',nfeatures_use,'_rd_sce_v2.rds'))
  
  plt_feature <- colnames(meta_info)[grepl('treat',colnames(meta_info))] #one of treat, or treatmentSt
  if(length(plt_feature)>0){
    p21 <- Seurat::DimPlot(srt, reduction = "umap", group.by = plt_feature)
    # p21 <- Seurat::DimPlot(srt, reduction = "pca", group.by = plt_feature)
    png(paste0(output_dir,"UMAP_",datatag,".png"), height = 2*400, width=2*600,res = 2*72)
    print(p21)
    dev.off()
  }
  
  
  # Seurat clusters results
  p22 <- Seurat::DimPlot(srt, reduction = "umap")
  # p21 <- Seurat::DimPlot(srt, reduction = "pca", group.by = plt_feature)
  png(paste0(output_dir,"UMAP_clusters_",datatag,"_12clusters.png"), height = 2*400, width=2*600,res = 2*72)
  print(p22)
  dev.off()
  
  # p <- cowplot::plot_grid(p21, p22)
  
  p23 <- Seurat::DimPlot(srt, reduction = "pca")
  # p21 <- Seurat::DimPlot(srt, reduction = "pca", group.by = plt_feature)
  png(paste0(output_dir,"PCA_clusters_",datatag,".png"), height = 2*400, width=2*600,res = 2*72)
  print(p23)
  dev.off()
  
  plt_clone <- colnames(meta_info)[grepl('clone',colnames(meta_info))] #one of treat, or treatmentSt
  if(length(plt_clone)>0){
    p24 <- Seurat::DimPlot(srt, reduction = "umap", group.by = plt_clone)
    # p21 <- Seurat::DimPlot(srt, reduction = "pca", group.by = plt_feature)
    png(paste0(output_dir,"UMAP_",plt_clone,'_',datatag,".png"), height = 2*400, width=2*400*length(plt_clone),res = 2*72)
    print(p24)
    dev.off()
  }
  top_hvg <- HVFInfo(srt) %>%
    mutate(., bc = rownames(.)) %>%
    arrange(desc(variance.standardized)) %>%
    top_n(nfeatures_use, variance.standardized)
  
  print(sum(top_hvg$bc %in% var_genes))
  
  data.table::fwrite(top_hvg, paste0(output_dir, datatag,'_',nfeatures_use, "_hvg_genes.csv"))

}  


plot_legend_treatmentSt <- function(save_dir){
  # cyan, corn flower, cobalt, dark blue, midnight blue
  rx_cols <- c('#00FFFF','#7aa4f0','#3a76c9', '#2a2aa3','#20208a')
  names(rx_cols) <- c('1Rx','2Rx','3Rx','4Rx','5Rx')
  #light grey, darkgray, dimgray, dark grey, dark dark grey
  # unrx_cols <- c('#989898','#676767','#696969', '#505050','#383838')
  unrx_cols <- c('#989898','#676767','#4d4d4d', '#2e2c2c','#0f0f0f')
  names(unrx_cols) <- c('1UnRx','2UnRx','3UnRx','4UnRx','5UnRx')
  
  rxh_cols <- c('#D5FF00','#FFFF00','#FFD500', '#FFBF00','#767600')
  names(rxh_cols) <- c('1RxH','2RxH','3RxH','4RxH','5RxH')
  
  cols <- c(unrx_cols, rx_cols, rxh_cols)
  lg <- plot_colors(cols, save_dir, legend_label='Treatment Cycles', 5)
  # lg
  return(lg)
}
set_color_treatmentSt_simplify_version <- function(clusters, datatag, cols=NULL){
  library(stringr)
  clusters_labels <- unique(clusters)
  metacells <- data.frame(treatmentSt=clusters_labels, stringsAsFactors=F)
  if(is.null(cols)){
    # cyan, corn flower, cobalt, dark blue, midnight blue
    rx_cols <- c('#00FFFF','#7aa4f0','#7aa4f0', '#7aa4f0','#7aa4f0') ## only keep first time point with contrast color
    names(rx_cols) <- c('1Rx','2Rx','3Rx','4Rx','5Rx')
    #light grey, darkgray, dimgray, dark grey, dark dark grey
    # unrx_cols <- c('#989898','#676767','#696969', '#505050','#383838')
    unrx_cols <- c('#989898','#4d4d4d','#4d4d4d', '#4d4d4d','#4d4d4d')
    names(unrx_cols) <- c('1UnRx','2UnRx','3UnRx','4UnRx','5UnRx')
    
    rxh_cols <- c('#D5FF00','#FFFF00','#FFFF00', '#FFFF00','#FFFF00')
    names(rxh_cols) <- c('1RxH','2RxH','3RxH','4RxH','5RxH')
    
    cols <- c(unrx_cols, rx_cols, rxh_cols)
  }
  order_labels <- data.frame(label=names(cols), 
                             color_code=cols,
                             order_st=seq(1:length(cols)))  # to fix the order of labels in plot
  
  
  metacells$label <- ifelse(grepl('T$',metacells$treatmentSt),'Rx','UnRx')
  metacells$label <- ifelse(grepl('TU$',metacells$treatmentSt),'RxH',metacells$label)
  if(datatag=='SA535'){
    metacells$label <- ifelse(metacells$label %in% c('Rx','RxH'),paste0(str_count(metacells$treatmentSt,'T'),
                                                                        metacells$label),
                              paste0(str_count(metacells$treatmentSt,'U')-1,metacells$label)) #SA535 contains 1 more U at starting points
  }else{
    metacells$label <- ifelse(metacells$label %in% c('Rx','RxH'),paste0(str_count(metacells$treatmentSt,'T'),
                                                                        metacells$label),
                              paste0(str_count(metacells$treatmentSt,'U'),metacells$label))
  }
  
  # cols1 <- cols[metacells$label]
  metacells <- metacells %>% dplyr::inner_join(order_labels, by='label')
  # names(cols1) <-  metacells$treatmentSt
  
  metacells <- metacells[order(metacells$order_st, decreasing=F),]
  # return(list(clone_palette=cols1,cols_use=cols1[clusters]))
  return(metacells)
  
}  

set_color_treatmentSt <- function(clusters, datatag, cols=NULL){
  library(stringr)
  clusters_labels <- unique(clusters)
  metacells <- data.frame(treatmentSt=clusters_labels, stringsAsFactors=F)
  if(is.null(cols)){
    # cyan, corn flower, cobalt, dark blue, midnight blue
    rx_cols <- c('#00FFFF','#7aa4f0','#3a76c9', '#2a2aa3','#20208a')
    names(rx_cols) <- c('1Rx','2Rx','3Rx','4Rx','5Rx')
    #light grey, darkgray, dimgray, dark grey, dark dark grey
    # unrx_cols <- c('#989898','#676767','#696969', '#505050','#383838')
    unrx_cols <- c('#989898','#676767','#4d4d4d', '#2e2c2c','#0f0f0f')
    names(unrx_cols) <- c('1UnRx','2UnRx','3UnRx','4UnRx','5UnRx')
    
    rxh_cols <- c('#D5FF00','#FFFF00','#FFD500', '#FFBF00','#767600')
    names(rxh_cols) <- c('1RxH','2RxH','3RxH','4RxH','5RxH')
    
    cols <- c(unrx_cols, rx_cols, rxh_cols)
  }
  order_labels <- data.frame(label=names(cols), 
                             color_code=cols,
                             order_st=seq(1:length(cols)))  # to fix the order of labels in plot
  
  
  metacells$label <- ifelse(grepl('T$',metacells$treatmentSt),'Rx','UnRx')
  metacells$label <- ifelse(grepl('TU$',metacells$treatmentSt),'RxH',metacells$label)
  if(datatag=='SA535'){
    metacells$label <- ifelse(metacells$label %in% c('Rx','RxH'),paste0(str_count(metacells$treatmentSt,'T'),
                                                                        metacells$label),
                              paste0(str_count(metacells$treatmentSt,'U')-1,metacells$label)) #SA535 contains 1 more U at starting points
  }else{
    metacells$label <- ifelse(metacells$label %in% c('Rx','RxH'),paste0(str_count(metacells$treatmentSt,'T'),
                                                                        metacells$label),
                              paste0(str_count(metacells$treatmentSt,'U'),metacells$label))
  }
  
  # cols1 <- cols[metacells$label]
  metacells <- metacells %>% dplyr::inner_join(order_labels, by='label')
  # names(cols1) <-  metacells$treatmentSt
  
  metacells <- metacells[order(metacells$order_st, decreasing=F),]
  # return(list(clone_palette=cols1,cols_use=cols1[clusters]))
  return(metacells)
  
}  

#' Assign a color to each cell based on some value
#' 
#' @param cell_vars Vector indicating the value of a variable associated with cells.
#' @param pal_fun Palette function that returns a vector of hex colors, whose
#' argument is the length of such a vector.
#' @param ... Extra arguments for pal_fun.
#' @return A vector of hex colors with one entry for each cell.
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}


# Based on monocle3 genes module method
get_genes_modules <- function(preprocess_mat, reduction_method = c("UMAP"), max_components = 10, 
                              umap.metric = "cosine", umap.min_dist = 0.1, umap.n_neighbors = 15L, 
                              umap.fast_sgd = FALSE, umap.nn_method = "annoy", k = 20, 
                              leiden_iter = 1, partition_qval = 0.05, weight = FALSE, 
                              resolution = NULL, random_seed = 0L, cores = 5, verbose = F, 
                              ...) {
  method = "leiden"
  print('Run reduction features...')
  if (random_seed != 0L) 
    set.seed(random_seed)
  
  
  umap_res = uwot::umap(as.matrix(preprocess_mat), n_components = max_components, 
                        metric = umap.metric, min_dist = umap.min_dist, n_neighbors = umap.n_neighbors, 
                        fast_sgd = umap.fast_sgd, n_threads = cores, verbose = verbose, 
                        nn_method = umap.nn_method, ...)
  row.names(umap_res) <- row.names(preprocess_mat)
  colnames(umap_res) <- paste0("dim_", 1:ncol(umap_res))
  reduced_dim_res <- umap_res
  if (verbose) 
    message("Running leiden clustering algorithm ...")
  
  print('Running leiden clustering algorithm ...')
  # print(row.names(t(umap_res))[1:5])
  cluster_result <- monocle3:::leiden_clustering(data = reduced_dim_res, 
                                                 pd = preprocess_mat, 
                                                 k = k, weight = weight, num_iter = leiden_iter, resolution_parameter = resolution, 
                                                 random_seed = random_seed, verbose = verbose, ...)
  # print('Running partitions algorithm ...')
  # cluster_graph_res <- monocle3:::compute_partitions(cluster_result$g, 
  #                                                    cluster_result$optim_res, partition_qval, verbose)
  # partitions <- igraph::components(cluster_graph_res$cluster_g)$membership[cluster_result$optim_res$membership]
  # names(partitions) <- row.names(reduced_dim_res)
  # partitions <- as.factor(partitions)
  print('Get genes modules ...')
  gene_module_df <- tibble::tibble(id = row.names(preprocess_mat), 
                                   module = factor(igraph::membership(cluster_result$optim_res))
                                   )#supermodule = partitions
  gene_module_df <- tibble::as_tibble(cbind(gene_module_df, 
                                            umap_res))
  return(gene_module_df)
  
  
}

get_treatment_passage_detail_desc <- function(metacell,datatag=NULL){
  # ts <- unique(treatmentSts)
  library(stringr)
  
  metacell$label <- ifelse(grepl('T$',metacell$treatmentSt),'Rx','UnRx')
  metacell$label <- ifelse(grepl('TU$',metacell$treatmentSt),'RxH',metacell$label)
  # unique(metacell$timepoint)
  metacell$treatment_desc <- paste0(metacell$label,'-', metacell$timepoint)
  # unique(metacell$treatment_desc)
  return(metacell)
}

get_treatment_detail_desc <- function(treatmentSts,datatag='SA609'){
  # ts <- unique(treatmentSts)
  library(stringr)
  metacells <- data.frame(treatmentSt=treatmentSts)
  metacells$label <- ifelse(grepl('T$',metacells$treatmentSt),'Rx','UnRx')
  metacells$label <- ifelse(grepl('TU$',metacells$treatmentSt),'RxH',metacells$label)
  if(datatag=='SA535'){
    metacells$label <- ifelse(metacells$label %in% c('Rx','RxH'),paste0(str_count(metacells$treatmentSt,'T'),
                                                                        metacells$label),
                              paste0(str_count(metacells$treatmentSt,'U')-1,metacells$label)) #SA535 contains 1 more U at starting points
  }else{
    metacells$label <- ifelse(metacells$label %in% c('Rx','RxH'),paste0(str_count(metacells$treatmentSt,'T'),
                                                                        metacells$label),
                              paste0(str_count(metacells$treatmentSt,'U'),metacells$label))
  }
  return(metacells$label)
}
get_treatment_desc <- function(treatmentSts){
  ts <- unique(treatmentSts)
  df <- tibble::tibble()
  for(obs_treatment in c('UnRx','Rx','RxH')){
    if(obs_treatment=='Rx'){
      conds <- grep('*T$', ts, value=T)
    } else if(obs_treatment=='RxH'){
      conds <- grep('*TU$', ts, value=T)
    }else{
      conds <- c('U',grep('*UU$', ts, value=T))
    }
    tmp <- data.frame(treatmentSt=conds, treatment_desc=rep(obs_treatment, length(conds)))
    df <- dplyr::bind_rows(df, tmp)
  }
  tmp <- data.frame(treatmentSt=treatmentSts)
  tmp <- tmp %>% inner_join(df,by='treatmentSt')
  print(dim(tmp)[1])
  print(length(treatmentSts))
  return(tmp$treatment_desc)
}
get_unique_clone_id <- function(clone_labels){
  set.seed(42) 
  cls <- sapply(strsplit(clone_labels,'_'), function(x){
    if(length(x)==1){
      return(x[1])
    }else if(length(x)==2){
      idx <- sample(c(1,2),1)
      return(x[idx])
    }else{
      idx <- sample(c(1,2,3),1)
      return(x[idx])
    }
    
  })
  return(as.character(cls))
}

plot_trajectory_clones_prevalence_SA1035 <- function(metacell, meta_lineage, output_dir, save_dir, datatag){
  save_figs_dir <- paste0(output_dir, 'figs_v3/')
  if(!dir.exists(save_figs_dir)){
    dir.create(save_figs_dir)
  }
  # sce$clone_v1 <- sce$clone
  # sce$clone <- ''
  # rownames(clone_df) <- clone_df$cell_id
  # cells_use <- intersect(clone_df$cell_id, colnames(sce))
  # sce[,cells_use]$clone <- clone_df[cells_use,'clone']
  # sce$clone <- ifelse(sce$clone=='','unassigned',sce$clone)
  
  # sum(clone_df$cell_id %in% metacell$cell_id)
  # metacell <- colData(sce) %>% as.data.frame()
  # dim(metacell)
  # unique(metacell$clone)
  # summary(as.factor(sce$clone))
  sce$clone <- get_unique_clone_id(sce$clone)
  
  # metacell$clone <- get_unique_clone_id(metacell$clone)
  # metacell$clone <- ifelse(metacell$clone=='unassigned','None',metacell$clone)
  unique(metacell$clone)
  unique(metacell$treatmentSt)
  # metacell$treatment_desc <- get_treatment_desc(metacell$treatmentSt)
  # metacell$treatment_desc <- get_treatment_detail_desc(metacell$treatmentSt)
  metacell <- get_treatment_passage_detail_desc(metacell)
  # unique(metacell$cluster_label)
  # table(metacell$treatment_desc,metacell$clone)
  metadata <- metacell %>%
    # dplyr::filter(clone!='unassigned')%>% # & cluster_label!=0
    dplyr::select(treatment_desc,clone,cell_id,treatmentSt,label)%>%
    # dplyr::mutate(treatment_clone=paste0(treatment_desc,' clone ',clone))%>%
    dplyr::mutate(treatment_clone=paste0(label,'-clone ',clone))%>%
    dplyr::mutate(clone=paste0('clone ',clone))
  
  pseudo <- data.table::fread(paste0(save_dir,'slingshot_output/slingshot_SA1035_0_PCA_pseudotime.csv')) %>% as.data.frame()
  weight <- data.table::fread(paste0(save_dir,'slingshot_output/slingshot_SA1035_0_PCA_cellWeights.csv')) %>% as.data.frame()
  
  # dim(pseudo)
  # head(pseudo)
  weight <- weight %>% inner_join(metadata,by='cell_id')
  dim(weight)  
  unique(weight$treatment_clone)
  colnames(weight)
  # SA1035
  # meta_lineage <- data.frame(lineage=c("Lineage1","Lineage2","Lineage3","Lineage4"),
  #                            lineage_desc=c("Lineage2","Lineage4",
  #                                           "Lineage1","Lineage3"))  # version 1
  # meta_lineage <- data.frame(lineage=c("Lineage1","Lineage2","Lineage3","Lineage4","Lineage5"),
  #                            lineage_desc=c("L4: Rx,RxH","L3: UnRx, RxH",
  #                                           "L2: Rx, RxH","L1: 1Rx, Rx", "L5:Rx, UnRx"))
  rownames(meta_lineage) <- meta_lineage$lineage
  # plttlts <- c("L2: Rx,RxH","L4: RxH, UnRx",
  #              "L1: 1Rx,Rx","L3: Rx")
  # weight$treatment_desc[1:5]
  lgs <- grep('Lineage',colnames(weight),value=T)
  lgs <- lgs[lgs!='Lineage5'] # very small lineage, not important, remove from the considered list of lineages
  stat <- tibble::tibble()
  for(lg in lgs){
    ps <- weight %>%
      dplyr::filter(!!sym(lg)>0)
    
    ps <- ps %>%
      group_by(treatment_desc) %>%
      summarise(fraction_cells=n()/dim(ps)[1])%>%
      filter(fraction_cells>0.02)
    # ps <- ps %>%
    #   group_by(treatment_clone) %>%
    #   summarise(fraction_cells=n()/dim(ps)[1])#%>%
    #filter(fraction_cells>0.1)
    ps$lineage <- lg
    stat <- dplyr::bind_rows(stat, ps)
  }
  stat1 <- stat %>%
    tidyr::pivot_wider(names_from='treatment_desc', 
                       values_from='fraction_cells',values_fill = 0)%>%
    tibble::column_to_rownames('lineage')%>%
    t()  
  # rns <- c("1Rx","2Rx","3Rx","4Rx",
  #          "1RxH","2RxH","3RxH",
  #          "1UnRx","2UnRx","3UnRx","4UnRx","5UnRx")
  # stat1 <- stat1[rns,]
  # treatment_order <- data.table::fread(paste0(save_figs_dir,'treatment_meta_colors.csv')) %>% as.data.frame()
  # treatment_order <- treatment_order %>%
  #   dplyr::filter(label %in% rownames(stat1)) %>%
  #   dplyr::arrange(!desc(order_st))
  # if(dim(treatment_order)[1]!=dim(stat1)[1]){
  #   stop('Double check meta treatment table')
  # }
  # stat1 <- stat1[treatment_order$label,]
  # rownames(stat1)
  colnames(stat1) <- meta_lineage[colnames(stat1),'lineage_desc']
  # stat1[stat1<0.05] <- 0
  # max(stat1)
  # viz_hm_proportions(stat1, datatag, output_dir, 'treatment_desc',ht=400, wd=200, 0.1)
  pts <- viz_hm_proportions(stat1, datatag, save_figs_dir, tag='treatment_desc', 
                            ht=400, wd=200, verbose_thres=0.1, col_fun=NULL, text_size = 11)
  
  stat <- tibble::tibble()
  for(lg in lgs){
    ps <- weight %>%
      dplyr::filter(!!sym(lg)>0 & !grepl('unassigned',treatment_clone))
    # head(ps)
    # dim(weight)
    # dim(ps)
    # summary(ps$Lineage4)
    # ps <- ps %>%
    #   group_by(clone) %>%
    #   summarise(fraction_cells=n()/dim(ps)[1])
    # ps1 <- ps %>%
    #   group_by(treatment_clone) %>%
    #   summarise(fraction_cells=n())
    # # sum(ps1$fraction_cells)
    ps <- ps %>%
      group_by(treatment_clone) %>%
      summarise(fraction_cells=n()/dim(ps)[1])%>%
      filter(fraction_cells>=0.03)
    ps$lineage <- lg
    stat <- dplyr::bind_rows(stat, ps)
  }
  dim(stat)
  # stat$treatment_clone <- stat$clone
  # unique(stat$treatment_clone)
  # stat$treatment_clone <- ifelse(grepl('cloneH',stat$treatment_clone) 
  #                                & !grepl('UnRx: cloneH',stat$treatment_clone),
  #                                paste0(stat$treatment_clone,'(Res)'),stat$treatment_clone)
  # 
  # stat$treatment_clone <- ifelse(grepl('UnRx: cloneE',stat$treatment_clone),
  #                                paste0(stat$treatment_clone,'(Sen)'),stat$treatment_clone)
  # 
  # stat$treatment_clone <- ifelse(grepl('1RxH: cloneB',stat$treatment_clone) | 
  #                                grepl('2RxH: cloneG',stat$treatment_clone) | 
  #                                grepl('3RxH: cloneG',stat$treatment_clone)  ,
  #                                paste0(stat$treatment_clone,'(DH, Res)'),stat$treatment_clone)
  
  # 
  # stat$treatment_clone <- ifelse(grepl('clone H',stat$treatment_clone),paste0(stat$treatment_clone,' (resistant*)'),
  #                                stat$treatment_clone)
  # stat$treatment_clone <- ifelse(grepl('clone E',stat$treatment_clone),paste0(stat$treatment_clone,' (sensitive*)'),
  #                                stat$treatment_clone)
  
  
  # stat$treatment_clone <- ifelse(stat$treatment_clone=='RxH: cloneB','RxH: cloneB (reversibility*)',stat$treatment_clone)
  # 
  # stat$treatment_clone <- ifelse(grepl('cloneA',stat$treatment_clone),paste0(stat$treatment_clone,' (resistant*)'),
  #                                                                            stat$treatment_clone)
  # stat$treatment_clone <- ifelse(grepl('cloneH',stat$treatment_clone),paste0(stat$treatment_clone,' (sensitive*)'),
  #                                stat$treatment_clone)
  # unique(stat$treatment_clone)
  # stat$treatment_clone <- paste0(stat$lineage,': ',stat$treatment_clone)
  # View(stat)
  # View(stat1)
  stat1 <- stat %>%
    # dplyr::select(-clone)%>%
    tidyr::pivot_wider(names_from='treatment_clone', 
                       values_from='fraction_cells',values_fill = 0)%>%
    tibble::column_to_rownames('lineage')%>%
    t()
  
  # stat1 <- stat %>%
  #   # dplyr::select(-clone)%>%
  #   tidyr::pivot_wider(names_from='clone', 
  #                      values_from='fraction_cells',values_fill = 0)%>%
  #   tibble::column_to_rownames('lineage')%>%
  #   t()
  
  
  # SA1035: drug holiday:
  #   treatment: H, untreated: E,
  # drug holiday:X6: B, X7, X8: G
  rownames(stat1)
  colnames(stat1)
  
  colnames(stat1) <- meta_lineage[colnames(stat1),'lineage_desc']
  # col_fun = colorRamp2(c(0, 0.05, 0.07, 0.1, max(stat1)), 
  #                      c("white","#ccccff",'#6666ff', "#0000ff","#000099"))
  # viz_hm_proportions(stat1, datatag, output_dir, 'treatment_clone',ht=400, wd=250, 0.05,col_fun)
  vals <- c(0, 0.05, 0.07, 0.1, 0.15, max(stat1))
  col_fun_desc = colorRamp2(vals, 
                            colorRampPalette(c("white", "darkgreen"))(length(vals)))
  pc <- viz_hm_proportions(stat1, datatag, save_figs_dir, tag='treatment_clone', 
                           ht=400, wd=200, verbose_thres=0.04, col_fun=col_fun_desc, text_size = 7)
  return(list(pts=pts, pc=pc))
}


viz_hm_proportions <- function(stat1, datatag, output_dir, tag='', 
                               ht=400, wd=200, verbose_thres=0.05, 
                               col_fun=NULL, text_size=11){
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  # library(circlize)
  # library(ComplexHeatmap)
  stat1 <- as.data.frame(stat1)
  if(is.null(col_fun)){
    # col_fun = colorRamp2(c(0, 0.05, 0.1, 0.2, max(stat1)), 
    #                      c("white","#ccccff",'#6666ff', "#0000ff","#000099"))
    vals <- c(0, 0.05, 0.1, 0.15, 0.2, max(stat1))
    col_fun = colorRamp2(vals, 
                         colorRampPalette(c("white", "darkgreen"))(length(vals)))
  }
  
  cell_funcs = function(j, i, x, y, width, height, fill) {
    if(stat1[i, j] >= verbose_thres){
      grid.text(sprintf("%.2f", round(stat1[i, j],2)), x, y, gp = gpar(fontsize = text_size))  #round(stat1[i, j], 3)
    }
  }
  tmp <- data.frame(lbs=rownames(stat1))
  tmp$ts_status <- ifelse(grepl('UnRx',tmp$lbs),'UnRx',
                          ifelse(grepl('RxH',tmp$lbs),'RxH','Rx'))
  tmp <- tmp %>%
    arrange(desc(ts_status))
  stat1 <- stat1[tmp$lbs,]
  p <- ComplexHeatmap::Heatmap(as.matrix(stat1), na_col = "black",
                               show_column_names=T,
                               show_row_names = T,
                               cluster_rows=F,cluster_columns=F,
                               name = "Prevalence", 
                               col = col_fun,
                               # row_order = tmp$lbs,
                               # row_split= obs_genes_df,
                               # row_title_rot = 0,
                               row_gap = unit(2, "mm"),
                               # column_split = obs_cells_df, 
                               # column_title = "ddd",
                               column_gap = unit(2, "mm"),
                               column_names_gp = grid::gpar(fontsize = 11),
                               row_names_gp = grid::gpar(fontsize = 11),
                               show_heatmap_legend = T,
                               # top_annotation=anno_col,
                               # left_annotation = left_anno,
                               cell_fun = cell_funcs,
                               row_dend_reorder=F,
                               heatmap_legend_param = list(direction = "vertical")
                               # column_title = "Pseudotime", 
                               # column_title_side = "bottom"
  )
  # p
  
  # add title
  
  png(paste0(output_dir,datatag, "_",tag,"_trajectory_prevalence.png"), height = 2*ht, width=2*wd, res = 2*72)
  ComplexHeatmap::draw(p, annotation_legend_side = "bottom",
       heatmap_legend_side = "bottom",
       padding = unit(c(1, 1, 2, 2), "mm"))
  dev.off()
  saveRDS(p, paste0(output_dir,datatag, "_",tag,"_trajectory_prevalence.rds"))
  # png(paste0(output_dir,datatag, "trajectory_prevalence_treatment.png"), height = 2*400, width=2*300, res = 2*72)
  # draw(p, annotation_legend_side = "bottom",
  #      heatmap_legend_side = "bottom",
  #      padding = unit(c(1, 1, 1, 1), "mm"))
  # dev.off()
  return(p)
}


plot_trajectory_clones_prevalence_SA609 <- function(metacell, input_dir, output_dir, meta_lineage=NULL){
  # For SA609
  # G, H
  # B, C
  # R
  # D
  if(is.null(meta_lineage)){
    # meta_lineage <- data.frame(lineage=c("Lineage1","Lineage2","Lineage3"),
    #                            lineage_desc=c("Lineage 3: Rx","Lineage 2: RxH","Lineage 1: UnRx"))
    meta_lineage <- data.frame(lineage=c("Lineage1","Lineage2","Lineage3"),
                               lineage_desc=c("L3-Rx","L2-RxH","L1-UnRx"))
  }
  rownames(meta_lineage) <- meta_lineage$lineage  
  
  
  # dim(metacell)
  # unique(metacell$clone)
  metacell$clone <- get_unique_clone_id(metacell$clone)
  # unique(metacell$treatmentSt)
  
  ## SA609, new labels, replace R by A 
  metacell <- metacell %>%
    dplyr::mutate(clone=replace(clone,clone=='R','A'))
  # unique(metacell$treatmentSt)
  # metacell$treatment_desc <- get_treatment_desc(metacell$treatmentSt)
  # metacell$treatment_desc <- get_treatment_detail_desc(metacell$treatmentSt)
  metacell <- get_treatment_passage_detail_desc(metacell, datatag='SA609')
  
  metadata <- metacell %>%
    dplyr::select(treatment_desc,clone,cell_id,treatmentSt,label)%>%
    # dplyr::mutate(treatment_clone=paste0(treatment_desc,' clone ',clone))%>%
    dplyr::mutate(treatment_clone=paste0(label,'-clone ',clone))%>%
    dplyr::mutate(clone=paste0('clone ',clone))
  

  ## Loading slingshot output from file
  pseudo <- data.table::fread(paste0(save_dir,'slingshot_SA609_10_PCA_pseudotime.csv')) %>% as.data.frame()
  weight <- data.table::fread(paste0(save_dir,'slingshot_SA609_10_PCA_cellWeights.csv')) %>% as.data.frame()
  
  ## Loading slingshot output directly from objects
  # pseudo <- slingPseudotime(crv_umap_embed) %>% as.data.frame()
  # weight <- slingshot::slingCurveWeights(crv_umap_embed) %>% as.data.frame()
  # dim(weight)
  # weight$cell_id <- rownames(weight)
  dim(weight)  
  sum(weight$Lineage1==0)
  pseudo$cell_id[1]
  weight <- weight %>% inner_join(metadata,by='cell_id')
  dim(weight)  
  # weight$treatment_desc[1:5]
  lgs <- grep('Lineage',colnames(weight),value=T)
  
  # Plot treatment desc heatmap
  stat <- tibble::tibble()
  for(lg in lgs){
    ps <- weight %>%
      dplyr::filter(!!sym(lg)>0)
    
    ps <- ps %>%
      group_by(treatment_desc) %>%
      summarise(fraction_cells=n()/dim(ps)[1])
    # ps <- ps %>%
    #   group_by(treatment_clone) %>%
    #   summarise(fraction_cells=n()/dim(ps)[1])#%>%
    #filter(fraction_cells>0.1)
    ps$lineage <- lg
    stat <- dplyr::bind_rows(stat, ps)
  }
  stat1 <- stat %>%
    tidyr::pivot_wider(names_from='treatment_desc', 
                       values_from='fraction_cells',values_fill = 0)%>%
    tibble::column_to_rownames('lineage')%>%
    t()  
  
  # rns <- c("1Rx","2Rx","3Rx","4Rx","1RxH","2RxH","3RxH","1UnRx","2UnRx","3UnRx","4UnRx","5UnRx")
  # treatment_order <- data.table::fread(paste0(output_dir,'treatment_meta_colors.csv')) %>% as.data.frame()
  # treatment_order <- treatment_order %>%
  #   dplyr::filter(label %in% rownames(stat1)) %>%
  #   dplyr::arrange(-order_st)
  # if(dim(treatment_order)[1]!=dim(stat1)[1]){
  #   stop('Double check meta treatment table')
  # }
  
  # stat1 <- stat1[treatment_order$label,]
  # stat1[stat1<0.05] <- 0
  # View(stat1)
  # max(stat1)
  
  # plttlts <- c("Lineage 3: Rx,RxH","Lineage 2: Rx,RxH",
  #              "Lineage 4: Rx,1Rx","Lineage 1: UnRx")
  # stat1[stat1<0.05] <- 0
  # View(stat1)
  # max(stat1)
  # dim(stat1)
  
  colnames(stat1) <- meta_lineage[colnames(stat1),'lineage_desc']
  
  pts <- viz_hm_proportions(stat1, datatag, paste0(output_dir,'figs_v3/'), tag='treatment', 
                          ht=400, wd=200, verbose_thres=0.05, col_fun=NULL)
  # pts
  # Plot clone proportion heatmap
  stat <- tibble::tibble()
  for(lg in lgs){
    ps <- weight %>%
      dplyr::filter(!!sym(lg)>0)
    
    # ps <- ps %>%
    #   group_by(clone) %>%
    #   summarise(fraction_cells=n()/dim(ps)[1])%>%
    #   filter(fraction_cells>0.01)
    ps <- ps %>%
      group_by(treatment_clone) %>%
      summarise(fraction_cells=n()/dim(ps)[1]) %>%
      filter(fraction_cells>0.01)
    ps$lineage <- lg
    stat <- dplyr::bind_rows(stat, ps)
  }
  # dim(stat)
  # stat$treatment_clone <- stat$clone
  # stat$treatment_clone <- ifelse(stat$treatment_clone=='clone B','clone B (*Rx)',stat$treatment_clone)
  
  # stat$treatment_clone <- ifelse(grepl('clone A',stat$treatment_clone),
  #                                paste0(stat$treatment_clone,' (*Rx)'),
  #                                stat$treatment_clone)
  # stat$treatment_clone <- ifelse(grepl('clone H',stat$treatment_clone),
  #                                paste0(stat$treatment_clone,' (*UnRx)'),
  #                                stat$treatment_clone)
  # stat$treatment_clone <- ifelse(stat$treatment_clone=='RxH: cloneB','RxH: cloneB (reversibility*)',stat$treatment_clone)
  # 
  # stat$treatment_clone <- ifelse(grepl('cloneA',stat$treatment_clone),paste0(stat$treatment_clone,' (resistant*)'),
  #                                                                            stat$treatment_clone)
  # stat$treatment_clone <- ifelse(grepl('cloneH',stat$treatment_clone),paste0(stat$treatment_clone,' (sensitive*)'),
  #                                stat$treatment_clone)
  print(unique(stat$treatment_clone))
  # stat$treatment_clone <- paste0(stat$lineage,': ',stat$treatment_clone)
  # View(stat)
  # View(stat1)
  stat1 <- stat %>%
    # dplyr::select(-clone)%>%
    tidyr::pivot_wider(names_from='treatment_clone', 
                       values_from='fraction_cells',values_fill = 0)%>%
    tibble::column_to_rownames('lineage')%>%
    t()
  
  # max(stat1)
  # dim(stat1)
  # rownames(stat1)
  colnames(stat1) <- meta_lineage[colnames(stat1),'lineage_desc']
  p <- viz_hm_proportions(stat1, datatag, paste0(output_dir,'figs_v3/'), tag='clone', 
                          ht=400, wd=200, verbose_thres=0.1, col_fun=NULL)
  return(p)
}


plot_trajectory_clones_prevalence_SA535 <- function(metacell, input_dir, output_dir, meta_lineage=NULL){
  
  save_figs_dir <- paste0(output_dir,'figs_v3/')
  if(!dir.exists(save_figs_dir)){
    dir.create(save_figs_dir)
  }  
  if(is.null(meta_lineage)){
    # meta_lineage <- data.frame(lineage=c("Lineage1","Lineage2","Lineage3","Lineage4"),
    #                            lineage_desc=c("L3: UnRx,RxH","L2: Rx,RxH",
    #                                           "L4: Rx,1Rx","L1: UnRx"))
    meta_lineage <- data.frame(lineage=c("Lineage1","Lineage2","Lineage3","Lineage4"),
                               lineage_desc=c("L3-Rx,RxH","L2-Rx,RxH",
                                              "L4-1Rx","L1-UnRx"))
  }
  rownames(meta_lineage) <- meta_lineage$lineage
  # grep('clone',colnames(metacell), value=T)
  unique(metacell$treatmentSt)
  metacell$clone <- get_unique_clone_id(metacell$clone)
  metacell$treatmentSt[1]
  # metacell$treatment_desc <- get_treatment_desc(metacell$treatmentSt)
  # metacell$treatment_desc <- get_treatment_detail_desc(metacell$treatmentSt, datatag = 'SA535')
  metacell <- get_treatment_passage_detail_desc(metacell,datatag='SA535')
  # unique(metacell$treatment_desc)
  # table(metacell$treatment_desc,metacell$clone)
  # metacell$cell_id <- paste0(metacell$library_id,'_',metacell$Barcode)
  # metadata <- metacell %>%
  #   dplyr::select(treatment_desc,cell_id,treatmentSt)
  
  metadata <- metacell %>%
    # dplyr::filter(clone!='unassigned')%>%
    dplyr::select(treatment_desc,clone,cell_id,treatmentSt,label)%>%
    # dplyr::mutate(treatment_clone=paste0(treatment_desc,' clone ',clone))%>%
    dplyr::mutate(treatment_clone=paste0(label,'-clone ',clone))%>%
    dplyr::mutate(clone=paste0('clone ',clone))
  pseudo <- data.table::fread(paste0(save_dir,'slingshot_SA535_7_PCA_pseudotime.csv')) %>% as.data.frame()
  weight <- data.table::fread(paste0(save_dir,'slingshot_SA535_7_PCA_cellWeights.csv')) %>% as.data.frame()
  # dim(weight)  
  
  
  weight <- weight %>% inner_join(metadata, by='cell_id')
  # dim(weight)  
  lgs <- grep('Lineage',colnames(weight),value=T)
  stat <- tibble::tibble()
  for(lg in lgs){
    ps <- weight %>%
      dplyr::filter(!!sym(lg)>0)
    
    ps <- ps %>%
      group_by(treatment_desc) %>%
      summarise(fraction_cells=n()/dim(ps)[1])
    # ps <- ps %>%
    #   group_by(treatment_clone) %>%
    #   summarise(fraction_cells=n()/dim(ps)[1])#%>%
    #filter(fraction_cells>0.1)
    ps$lineage <- lg
    stat <- dplyr::bind_rows(stat, ps)
  }
  
  
  stat1 <- stat %>%
    tidyr::pivot_wider(names_from='treatment_desc', 
                       values_from='fraction_cells',values_fill = 0)%>%
    tibble::column_to_rownames('lineage')%>%
    t()  
  
  # treatment_order <- data.table::fread(paste0(save_figs_dir,'treatment_meta_colors.csv')) %>% as.data.frame()
  # treatment_order <- treatment_order %>%
  #   dplyr::filter(label %in% rownames(stat1)) %>%
  #   dplyr::arrange(-order_st)
  # if(dim(treatment_order)[1]!=dim(stat1)[1]){
  #   stop('Double check meta treatment table')
  # }
  # stat1 <- stat1[treatment_order$label,]
  # rownames(stat1)
  # rns <- c("1Rx","2Rx","3Rx","1RxH","2RxH","4RxH","1UnRx","2UnRx","3UnRx","4UnRx")
  # stat1 <- stat1[rns,]
  # colnames(stat1)
  # max(stat1)
  # dim(stat1)
  colnames(stat1) <- meta_lineage[colnames(stat1),'lineage_desc']
  # viz_hm_proportions(stat1, datatag, save_figs_dir, tag='treatment', ht=400, wd=190, col_fun=NULL)
  pts <- viz_hm_proportions(stat1, datatag, save_figs_dir, tag='treatment', 
                            ht=400, wd=200, verbose_thres=0.1, col_fun=NULL, text_size = 9)
  
  
  stat <- tibble::tibble()
  for(lg in lgs){
    ps <- weight %>%
      dplyr::filter(!!sym(lg)>0 & !grepl('unassigned',treatment_clone))
    
    # ps <- ps %>%
    #   group_by(clone) %>%
    #   summarise(fraction_cells=n()/dim(ps)[1])
    
    ps <- ps %>%
      group_by(treatment_clone) %>%
      summarise(fraction_cells=n()/dim(ps)[1])%>%
      filter(fraction_cells>0.02)
    ps$lineage <- lg
    stat <- dplyr::bind_rows(stat, ps)
  }
  # dim(stat)
  # SA535: treatment: A, untreated: G 
  # drug holiday:X8: D, X10: E
  # stat$treatment_clone <- stat$clone
  # unique(stat$treatment_clone)
  # stat$treatment_clone <- ifelse(stat$treatment_clone=='4RxH: cloneE','4RxH: cloneE(*Rx)',stat$treatment_clone)
  # stat$treatment_clone <- ifelse(stat$treatment_clone=='2RxH: cloneD','2RxH: cloneD(*Rx)',stat$treatment_clone)
  # stat$treatment_clone <- ifelse(grepl('cloneA',stat$treatment_clone),paste0(stat$treatment_clone,'(*Rx)'),
  #                                stat$treatment_clone)
  # stat$treatment_clone <- ifelse(grepl('UnRx: cloneG',stat$treatment_clone),paste0(stat$treatment_clone,'(*UnRx)'),
  #                                stat$treatment_clone)
  # stat$treatment_clone <- ifelse(stat$treatment_clone=='RxH: cloneB','RxH: cloneB (reversibility*)',stat$treatment_clone)
  # 
  # stat$treatment_clone <- ifelse(grepl('cloneA',stat$treatment_clone),paste0(stat$treatment_clone,' (resistant*)'),
  #                                                                            stat$treatment_clone)
  # stat$treatment_clone <- ifelse(grepl('cloneH',stat$treatment_clone),paste0(stat$treatment_clone,' (sensitive*)'),
  #                                stat$treatment_clone)
  # unique(stat$treatment_clone)
  # stat$treatment_clone <- paste0(stat$lineage,': ',stat$treatment_clone)
  # View(stat)
  # View(stat1)
  stat1 <- stat %>%
    # dplyr::select(-clone)%>%
    tidyr::pivot_wider(names_from='treatment_clone', 
                       values_from='fraction_cells',values_fill = 0)%>%
    tibble::column_to_rownames('lineage')%>%
    t()
  # max(stat1)
  vals <- c(0, 0.03, 0.05, 0.1, max(stat1))
  col_fun_desc = colorRamp2(vals, 
                       colorRampPalette(c("white", "darkgreen"))(length(vals)))
  
  # col_fun_desc = colorRamp2(c(0, 0.05, 0.07, 0.1, max(stat1)), 
  #                      c("white","#ccccff",'#6666ff', "#0000ff","#000099"))
  # viz_hm_proportions(stat1, datatag, output_dir, tag='treatment_clone', ht=400, wd=270, col_fun=col_fun_desc)
  colnames(stat1) <- meta_lineage[colnames(stat1),'lineage_desc']
  colnames(stat1) 
  rownames(stat1)[1]
  pc <- viz_hm_proportions(stat1, datatag, save_figs_dir, tag='treatment_clone', 
                          ht=400, wd=200, verbose_thres=0.03, col_fun=col_fun_desc, text_size = 7)
  # pc
  
    
}
# plot_lineages_oldversion <- function(sds, rd1, save_dir){
#   library(viridis)
#   nc <- 3
#   pt <- slingPseudotime(sds)
#   head(pt)
#   nms <- colnames(pt)
#   nr <- ceiling(length(nms)/nc)
#   pal <- viridis(100, end = 0.95)
#   png(paste0(save_dir,"slingshot_plot_lineages_pseudotime.png"), height = 2*800, width=2*900,res = 2*72)
#   par(mfrow = c(nr, nc))
#   for (i in nms) {
#     colors <- pal[cut(pt[,i], breaks = 100)]
#     plot(rd1, col = colors, pch = 16, cex = 0.5, main = i)
#     lines(SlingshotDataSet(sds), lwd = 1.5, col = 'black', type = 'lineages')
#   }
#   dev.off()
#   
#   t <- pt[,i]
#   class(pt)
# }


# plot_lineages_colorby_treatmentSt <- function(sds, rd1, meta_info, datatag, save_dir){
#   
#   nc <- 3
#   pt <- slingPseudotime(sds) %>% as.data.frame()
#   nms <- colnames(pt)
#   pt$cell_id <- as.character(rownames(pt))
#   # head(pt)
#   dim(pt)
#   class(pt)
#   meta_info$cell_id <- rownames(meta_info)
#   pt <- pt %>% inner_join(meta_info, by = c('cell_id'))
#   pt$treat
#   nr <- ceiling(length(nms)/nc)
#   res2 <- set_color_treatmentSt(pt[,'treatmentSt'], datatag, cols=NULL)
#   # pal <- viridis(100, end = 0.95)
#   png(paste0(save_dir,"slingshot_plot_lineages_pseudotime_treatmentSt.png"), height = 2*800, width=2*900,res = 2*72)
#   par(mfrow = c(nr, nc))
#   for (i in nms) {
#     # colors <- pal[cut(pt[,i], breaks = 100)]
#     tmp <- pt[,c(i,'treatmentSt')]
#     tmp$treatmentSt <- ifelse(!is.na(tmp[,i]),tmp$treat, NA)
#     colors <- res2$clone_palette[tmp$treatmentSt]
#     plot(rd1, col = colors, pch = 16, cex = 0.5, main = i)
#     lines(SlingshotDataSet(sds), lwd = 1.5, col = 'black', type = 'lineages')
#   }
#   dev.off()
# }  



# library(tidymodels)
# BiocManager::install("ranger")
# library("devtools")
# install_github("tidymodels/tidymodels")
# get_DE_genes_random_forests <- function(){
#   # Get top highly variable genes
#   top_hvg <- HVFInfo(seu) %>% 
#     mutate(., bc = rownames(.)) %>% 
#     arrange(desc(residual_variance)) %>% 
#     top_n(300, residual_variance) %>% 
#     pull(bc)
#   # Prepare data for random forest
#   dat_use <- t(logcounts(sce[top_hvg,]))
#   dat_use_df <- cbind(slingPseudotime(sds)[,4], dat_use) # Do curve 4, so 4st columnn
#   colnames(dat_use_df)[1] <- "pseudotime"
#   dat_use_df <- as.data.frame(dat_use_df[!is.na(dat_use_df[,1]),])
#   dat_split <- rsample::initial_split(dat_use_df)
#   dat_train <- rsample::training(dat_split)
#   dat_val <- rsample::testing(dat_split)
#   model <- rand_forest(mtry = 200, trees = 1400, min_n = 15, mode = "regression") %>%
#     set_engine("ranger", importance = "impurity", num.threads = 3) %>%
#     fit(pseudotime ~ ., data = dat_train)
#   
#   val_results <- dat_val %>% 
#     mutate(estimate = predict(model, .[,-1]) %>% pull()) %>% 
#     select(truth = pseudotime, estimate)
#   metrics(data = val_results, truth, estimate)
#   
#   var_imp <- sort(model$fit$variable.importance, decreasing = TRUE)
#   # View(head(var_imp))
#   top_genes <- names(var_imp)[1:6]
#   # Convert to gene symbol
#   # top_gene_name <- gns$gene_name[match(top_genes, gns$gene)]
#   genes_info <- data.table::fread("/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/Symbol_ensembl.csv") %>% as.data.frame() 
#   rownames(genes_info) <- genes_info$Ensembl
#   top_gene_name <- genes_info[top_genes,'Symbol']
#   
#   par(mfrow = c(3, 2))
#   for (i in seq_along(top_genes)) {
#     colors <- pal[cut(dat_use[,top_genes[i]], breaks = 100)]
#     plot(umaps, col = colors, 
#          pch = 16, cex = 0.5, main = top_gene_name[i])
#     lines(SlingshotDataSet(sds), lwd = 1.5, col = 'black', type = 'lineages')
#   }
#   return(var_imp)
#   
# }