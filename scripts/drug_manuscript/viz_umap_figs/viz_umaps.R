
suppressPackageStartupMessages({
  require("SingleCellExperiment")
  require("stringr")
  require("tidyverse")
  require("Seurat")
  # require("sctransform")
  require("dplyr")
  require("inlmisc")
})

library(extrafont)
font_import(prompt=F, paths ='/usr/share/fonts/truetype/myfonts/') # import Helvetica font
fonts()
my_font <- "Helvetica"
# Just modify function for predefined clone colors. 

make_clone_palette <- function(levels) {
  # install.packages("inlmisc", dependencies = TRUE)  # TO DO: check this package
  clone_names <- sort(levels)
  pal <- as.character(inlmisc::GetColors(length(clone_names)))  
  # if (length(levels) <= 12 & length(levels)>8) {
  #   pal <- brewer.pal(max(length(levels), 3), "Set3")
  # } else if (length(levels) <= 20 & length(levels) > 12) {
  #   pal <- clone_palette_20
  # } else {
  #   pal <- colorRampPalette(brewer.pal(8, "Set2"))(length(levels))
  #   print("WARNING: more clones than palette can accomodate!")
  # }
  pal[pal=='#E7EBFA']<- '#ABB9ED'
  names(pal) <- clone_names
  pal <- pal[levels]
  return(pal)
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
# basename <- 'SA1035'
# norm_sce <- readRDS('/home/htran/storage/datasets/drug_resistance/rna_results/rnaseq_v6/SA1035-v6/total_sce_v3.rds')
# dim(norm_sce)

# cell_clones <- read.csv('/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_total_v2/cell_clones.csv', check.names = F, stringsAsFactors = F)
# levels <- unique(cell_clones$clone_id)
get_color_clones <- function(tag, color_fn){
  # if(datatag=='SA609'){
  #   datatag <-'SA1000'
  # }
  # input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/'
  # color_df <- data.table::fread(paste0(input_dir,'all_clones_new_colors.csv')) %>% as.data.frame()
  # # color_df$cache_path <- NULL
  # color_df <- color_df %>%
  #   dplyr::filter(family == grep(paste0("*",datatag,"*"),unique(color_df$ref_datatag), value=T))#%>%
  # 
  # color_df <- color_df %>%
  #   dplyr::filter(family == grep(paste0("*",datatag,"*"),unique(color_df$ref_datatag), value=T))
  #   # dplyr::select(clone_id, datatag, family, colour)
  # unique(color_df$ref_datatag)
  
  color_df <- data.table::fread(color_fn)
  # color_df <- data.table::fread(paste0(output_dir,'colorcode_total.csv'))
  color_df <- color_df %>%
      dplyr::filter(datatag==tag)
  
  if(dim(color_df)[1] > 0){
    cols_use <- color_df$colour
    names(cols_use) <- color_df$clone_id
    cols_use['None'] <- '#D3D3D3'
    # if('None' %in% names(cols_use)){
    #     
    # }
    
  }else{
    # cols_use <- make_clone_palette(obs_clones)
    stop('Error, check color mapping')
  }
  return(cols_use)
}

thesis_theme <- ggplot2::theme(
  text = element_text(color="black",size = 8, hjust = 0.5, family=my_font),
  # axis.title.x = element_text(color="black",size=8, hjust = 0.5, family=my_font),
  # axis.title.y = element_text(color="black",size=8, hjust = 0.5, family=my_font),
  # axis.text.x = element_text(color="black",size=7, hjust = 0.5, family=my_font, angle = 90),
  # axis.text.x = element_blank(),
  axis.ticks = element_blank(),
  strip.placement = "outside",
  # axis.line = element_line(colour = "black"),
  axis.line = element_blank(),
  axis.text = element_blank(),
  axis.title = element_blank(),
  plot.title = element_text(color="black",size=11, face="bold",family=my_font, hjust = 0.5),
  # legend.title=element_text(color="black",size=7, hjust = 0.5, family=my_font),
  # legend.text=element_text(color="black",size=7, hjust = 0.5, family=my_font),
  strip.text.x = element_text(color="black",size=11, family=my_font),
  strip.text.y = element_text(color="black",size=11, family=my_font),
  # legend.position = lg_pos,
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(-2,-2,-2,-2),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "grey50", fill=NA, size=0.5),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank()
)

viz_umap_obs_clones <- function(umap_df, cols_use, datatag, output_dir, 
                                obs_treatment, obs_passage=NULL, plottitle='', plotlegend=F){
  xl <- c(min(umap_df$UMAP_1),max(umap_df$UMAP_1))
  yl <- c(min(umap_df$UMAP_2),max(umap_df$UMAP_2))
  
  my_font <- "Helvetica"
  if(!'treatmentSt' %in% colnames(umap_df)){
    if('treat' %in% colnames(umap_df)){
      umap_df <- umap_df %>%
        dplyr::rename(treatmentSt=treat)
    }else{
      stop('Please check input treatmentSt')
    }
  }
  
  umap_df <- umap_df %>%
    dplyr::mutate(clone=replace(clone, is.na(clone),'None'))
  # umap_df$treatmentSt <- umap_df$treat
  if(!is.null(obs_passage)){
    tmp <- umap_df %>%
      dplyr::filter(timepoint==obs_passage)
  }else{
    tmp <- umap_df
  }
  
  if(dim(tmp)[1]==0){
    print('Double check input passage info, exit')
    print('No cell at this passage')
    return(NULL)
  }
  ts <- unique(tmp$treatmentSt)
  if(obs_treatment=='Rx'){
    conds <- grep('*T$',ts, value=T)
    
  } else if(obs_treatment=='RxH'){
    conds <- grep('*TU$',ts, value=T)
  }else{
    conds <- c('U',grep('*UU$',ts, value=T))
  }
  tmp <- tmp %>%
    dplyr::filter(treatmentSt %in% conds)
  
  print(dim(tmp))
  print(obs_treatment)
  # print(unique(tmp$treatmentSt))
  
  # obs_clones=NULL
  obs_clones <- sort(unique(tmp$clone))
  unassign_clones <- c('unassigned','Unassigned','Un','un','None')
  obs_clones <- obs_clones[!obs_clones %in% unassign_clones]
  # print(obs_clones)
  # if(is.null(obs_clones)){
  #   
  # }
  cols_use <- cols_use[obs_clones]
  other_clones <- '#F0F0F0'
  names(other_clones) <- 'Others'
  cols_use <- c(cols_use,other_clones)
  # print(cols_use)
  # umap_df <- umap_df %>%
  #   dplyr::mutate(clone=replace(clone, !clone %in% obs_clones,'Others'))
  
  rownames(umap_df) <- umap_df$cell_id
  cells_excluded <- umap_df$cell_id[!umap_df$cell_id %in% tmp$cell_id]
  # umap_df[cells_excluded,'clone'] <- 'Others'
  # umap_df$alpha_val <- 1
  # umap_df[cells_excluded,'alpha_val'] <- 0.03
  # umap_df <- umap_df %>%
  #   dplyr::filter(cell_id %in% tmp$cell_id)
  # print(summary(as.factor(umap_df$clone)))
  
  
  umap_df$clone <- factor(umap_df$clone, levels = sort(unique(umap_df$clone)))
  p <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) + 
    geom_point(color='#e0e0e0') +  # grey color for all cells landscape displaying in background
    geom_point(data=tmp, aes(color=clone), size=0.7) + 
    scale_color_manual(values=cols_use, name='') + 
    labs(title = plottitle, x=' ') +
    xlim(xl[1], xl[2]) + 
    ylim(yl[1], yl[2])
    
  p <- p + thesis_theme
  # lg <- cowplot::get_legend(p + guides(color = guide_legend(nrow = 1, 
  #                                                           title.position = "left", 
  #                                                           override.aes = list(size=2))))
  # plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  if(!plotlegend){
    lg_pos <- "none"
    p <- p + ggplot2::theme(legend.position = lg_pos)  
  }else{
    lg_pos <- "right"
    p <- p + ggplot2::theme(legend.position = lg_pos)  
  }
  
  results <- list(p=p, cols_use=cols_use, df=tmp)
  # basename <- paste0(datatag,"_", gsub(' ','',plottitle))
  # saveRDS(results, paste0(output_dir, basename, ".rds"))
  # 
  # png(paste0(output_dir,basename,".png"), height = 2*350, width=2*500,res = 2*72)
  # print(p)
  # dev.off()
  return(results)
  
}

get_umap <- function(sce_fn, output_dir, datatag='SA'){
  sce <- readRDS(sce_fn)
  print(dim(sce))
  dims <-  1:30
  nb_hvg <- 3000
  # library(Seurat)
  # sce$library_id <- gsub('SCRNA10X_SA_CHIP','',sce$library_id)
  # sce$Grouping <- ifelse(sce$Site_origin=="Tumor_Recur","Primary",sce$Grouping)
  srt <- Seurat::as.Seurat(sce, counts = "counts", 
                           data="logcounts", assay = "RNA", project = "SingleCellExperiment") #set to NULL if only normalized data are present
  print(dim(srt))
  srt[["RNA"]]@scale.data <- as.matrix(logcounts(sce))
  srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = nb_hvg)
  # srt <- ScaleData(srt, features = all.genes) # scale for all genes, can scale by hvg genes
  srt <- RunPCA(object = srt, verbose = FALSE, npcs = 30)
  # srt <- JackStraw(srt, num.replicate = 100)
  # srt <- ScoreJackStraw(srt, dims = 1:20)
  # JackStrawPlot(srt, dims = 1:15)
  # ElbowPlot(srt)
  srt <- FindNeighbors(srt, dims = dims)
  # srt <- FindClusters(srt, resolution = res)
  srt <- RunUMAP(object = srt, dims = dims, verbose = FALSE) #umap.method = "umap-learn", metric='correlation'
  
  pca_df <- as.data.frame(Seurat::Embeddings(object = srt, reduction = "pca"))
  umap_df <- as.data.frame(Seurat::Embeddings(object = srt, reduction = "umap"))
  # library(dplyr)
  meta_info <- srt@meta.data
  meta_info$cell_id <- rownames(meta_info)
  # meta_info <- meta_info %>%
  #   dplyr::select(cell_id, Barcode, library_id, sample, treatmentSt, timepoint)
  
  umap_df$cell_id <- rownames(umap_df)
  umap_df <- umap_df %>% inner_join(meta_info, by=c("cell_id"))
  # umap_df$cell_id <- paste0(umap_df$library_id,'_',umap_df$Barcode)
  
  pca_df$cell_id <- rownames(pca_df)
  pca_df <- pca_df %>% inner_join(meta_info, by=c("cell_id"))
  # pca_df$cell_id <- paste0(pca_df$library_id,'_',pca_df$Barcode)
  print(dim(umap_df))
  print(dim(pca_df))
  
  data.table::fwrite(pca_df, paste0(output_dir, datatag, "_norm_pca.csv"))
  data.table::fwrite(umap_df, paste0(output_dir, datatag, "_norm_umap.csv"))
}


plot_fill_barplot_wholedataset_rnaseq <- function(umap_df, cols_use, 
                                              output_dir, datatag='SA', 
                                              plottitle=NULL, plotlegend=F){
  
  cns <- colnames(umap_df)
  tid <- cns[grepl('treatment',cns)]
  # if(nchar(tid)==0){
  #   tid <- cns[grepl('treat',cns)] ## sometimes column name is 'treat'
  #   if(is.null(tid)){
  #     print("BUG here, check treatment status")
  #   }  
  # }
  print(dim(umap_df))
  # sum(is.na(umap_df$clone_id))
  
  cid <- cns[grepl('clone',cns)]
  unassign_clones <- c('unassigned','Unassigned','Un','un','None')
  umap_df <- umap_df %>% 
    dplyr::rename(clone_id=!!sym(cid)) %>%
    dplyr::filter(!is.na(clone_id)) %>% 
    dplyr::mutate(clone_id = case_when(clone_id %in% unassign_clones ~ 'None',
                                      TRUE  ~ clone_id))%>% 
    dplyr::mutate(treatment_desc = case_when(grepl('TU$',!!sym(tid)) ~ 'RxH',
                                             grepl('T$',!!sym(tid)) ~ 'Rx',
                                             TRUE  ~ 'UnRx'))
  print(unique(umap_df$clone_id))
  print(unique(umap_df$treatment_desc))
  
  # fa='clone_id', xa='timepoint'
  cols_use <- cols_use[unique(umap_df$clone_id)]
  cls_df <- umap_df %>%
    dplyr::group_by(clone_id, timepoint, treatment_desc) %>%
    dplyr::summarise(Freq=n())%>%
    as.data.frame()
  # cls_df <- cls_df[gtools::mixedsort(cls_df$clone_id),]
  # dim(cls_df)
  
  
  cls_df$timepoint <- factor(cls_df$timepoint, levels=gtools::mixedsort(unique(cls_df$timepoint)))
  cls_df$clone_id <- factor(cls_df$clone_id, levels=sort(unique(cls_df$clone_id)))
  ts <- unique(cls_df$treatment_desc)
  ts1 <- ts[grepl('UnRx',ts)]
  ts2 <- ts[grepl('RxH',ts)]
  ts3 <- ts[!ts %in% c(ts1, ts2)]
  # cls_df$treatment_desc <- factor(cls_df$treatment_desc, levels=c(ts1, ts3, ts2))
  cls_df$treatment_desc <- paste0('10x-',cls_df$treatment_desc)
  cls_df$treatment_desc <- factor(cls_df$treatment_desc, levels=paste0('10x-',c(ts1, ts3, ts2)))
  my_font <- "Helvetica"
  p <- ggplot(cls_df, aes(fill=clone_id, y=Freq, x=timepoint)) + 
    geom_bar(position="fill", stat="identity",width=0.4) + 
    facet_grid(. ~ treatment_desc) + #, space='free', drop=F, scales="free".   . ~ treatment_desc
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) + 
    scale_fill_manual(values = cols_use)#name='Clone '
  p <- p + labs(title=plottitle, x=NULL, y=NULL) #x=xlabel, y=ylabel, y='Clonal fraction'
  
  thesis_theme <- ggplot2::theme(
    # text = element_text(color="black",size = 8, hjust = 0.5, family=my_font),
    # axis.title.x = element_text(color="black",size=8, hjust = 0.5, family=my_font),
    # axis.title.y = element_text(color="black",size=8, hjust = 0.5, family=my_font),
    axis.text.x = element_text(color="black",size=10, hjust = 0.5, family=my_font), #, angle = 90
    axis.text.y = element_text(color="black",size=8, hjust = 0.5, family=my_font), #, angle = 90
    # axis.text.x = element_blank(),
    # axis.ticks.y = element_blank(),
    strip.placement = "outside",
    # axis.line = element_line(colour = "black"),
    # axis.line.x = element_line(colour = "black"),
    # axis.line.y = element_line(colour = "black"),
    # axis.text.y = element_blank(),
    # axis.title = element_blank(),
    plot.title = element_text(color="black",size=11, face="bold",family=my_font, hjust = 0.5),
    legend.title=element_blank(),
    # legend.title=element_text(color="black",size=9, hjust = 0.5, family=my_font),
    legend.text=element_text(color="black",size=9, hjust = 0.5, family=my_font),
    legend.key.height= unit(0.2, 'cm'),
    legend.key.width= unit(0.5, 'cm'),
    strip.text.x = element_text(color="black",size=11, family=my_font, face="bold"),
    strip.text.y = element_text(color="black",size=11, family=my_font),
    strip.background = element_blank(),
    # legend.position = 'none',
    # panel.background = element_blank(),
    # panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )
  p <- p + theme_bw() +  thesis_theme
  lg <- cowplot::get_legend(p + guides(fill = guide_legend(nrow=3, title.position = "left", 
                                                           override.aes = list(shape = 0, size=0.2)))) #nrow=1
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  if(!plotlegend){
    lg_pos <- "none"
    p <- p + ggplot2::theme(legend.position = lg_pos)  
  }else{
    # lg_pos <- "right"
    lg_pos <- "top"
    p <- p + ggplot2::theme(legend.position = lg_pos)  
  }
  
  results <- list(plg=plg, p=p, cols_use=cols_use, df=cls_df, umap_df=umap_df)
  # saveRDS(results, paste0(output_dir, datatag, '_',plottitle,"_dlp_prevalence.rds"))
  # nrow <- length(unique(cls_df$timepoint))
  # png(paste0(output_dir,datatag, "_",plottitle, ".png"), height = 2*250, width=2*(50*nrow+20),res = 2*72)
  # print(p)
  # dev.off()
  # 
  # png(paste0(output_dir,datatag, "_",plottitle, "_legend.png"), height = 2*50, 
  #     width=2*(30*length(unique(cls_df$clone_id))+100),res = 2*72)
  # print(plg)
  # dev.off()
  
  return(results)
}

plot_clone_color_legend <- function(datatag, base_dir, ncols_grid=5){
  color_fn <- paste0(base_dir,'materials/umap_figs/colorcode_total_v3.csv.gz')
  cols_use <- get_color_clones(datatag, color_fn) # predefined clone color code for DLP and 10x inferred clones clonealign
  cols_use <- cols_use[gtools::mixedorder(cols_use)]
  cls_df <- data.frame(clone_id=names(cols_use))
  cls_df$Freq <- 1
  
  p <- ggplot(cls_df, aes(fill=clone_id, y=Freq, x=clone_id)) + 
    geom_bar(position="fill", stat="identity",width=0.4) + 
    scale_fill_manual(values = cols_use)+
    theme(legend.text=element_text(color="black",size=10, hjust = 0.5, family=my_font),
          legend.title=element_blank())
  lg <- cowplot::get_legend(p + guides(fill = guide_legend("Clone ",ncol=ncols_grid, title.position = "left", 
                                                           override.aes = list(shape = 0, size=0.1)))) #nrow=1
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  # plg
  return(plg)
}
plot_fill_barplot_wholedataset_v2 <- function(cell_clones, cols_use, 
                                           output_dir, datatag='SA', plottitle=NULL, plotlegend=F){
  
  # fa='clone_id', xa='timepoint'
  cols_use <- cols_use[unique(cell_clones$clone_id)]
  cls_df <- cell_clones %>%
    dplyr::group_by(clone_id, timepoint, treatment_desc) %>%
    dplyr::summarise(Freq=n())%>%
    as.data.frame()
  # cls_df <- cls_df[gtools::mixedsort(cls_df$clone_id),]
  # dim(cls_df)
  
  cls_df$timepoint <- factor(cls_df$timepoint, levels=gtools::mixedsort(unique(cls_df$timepoint)))
  cls_df$clone_id <- factor(cls_df$clone_id, levels=sort(unique(cls_df$clone_id)))
  ts <- unique(cls_df$treatment_desc)
  ts1 <- ts[grepl('UnRx',ts)]
  ts2 <- ts[grepl('RxH',ts)]
  ts3 <- ts[!ts %in% c(ts1, ts2)]
  cls_df$treatment_desc <- factor(cls_df$treatment_desc, levels=c(ts1, ts3, ts2))
  my_font <- "Helvetica"
  p <- ggplot(cls_df, aes(fill=clone_id, y=Freq, x=timepoint)) + 
    geom_bar(position="fill", stat="identity",width=0.4) + 
    facet_grid(. ~ treatment_desc) + #, space='free', drop=F, scales="free".   . ~ treatment_desc
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) + 
    scale_fill_manual(values = cols_use)#name='Clone '
  p <- p + labs(title=plottitle, x=NULL, y=NULL) #x=xlabel, y=ylabel, y='Clonal fraction'
  
  thesis_theme <- ggplot2::theme(
    # text = element_text(color="black",size = 8, hjust = 0.5, family=my_font),
    # axis.title.x = element_text(color="black",size=8, hjust = 0.5, family=my_font),
    # axis.title.y = element_text(color="black",size=8, hjust = 0.5, family=my_font),
    axis.text.x = element_text(color="black",size=10, hjust = 0.5, family=my_font), #, angle = 90
    axis.text.y = element_text(color="black",size=8, hjust = 0.5, family=my_font), #, angle = 90
    # axis.text.x = element_blank(),
    # axis.ticks.y = element_blank(),
    strip.placement = "outside",
    # axis.line = element_line(colour = "black"),
    # axis.line.x = element_line(colour = "black"),
    # axis.line.y = element_line(colour = "black"),
    # axis.text.y = element_blank(),
    # axis.title = element_blank(),
    plot.title = element_text(color="black",size=11, face="bold",family=my_font, hjust = 0.5),
    legend.title=element_blank(),
    # legend.title=element_text(color="black",size=9, hjust = 0.5, family=my_font),
    legend.text=element_text(color="black",size=9, hjust = 0.5, family=my_font),
    legend.key.height= unit(0.2, 'cm'),
    legend.key.width= unit(0.5, 'cm'),
    strip.text.x = element_text(color="black",size=11, family=my_font, face="bold"),
    strip.text.y = element_text(color="black",size=11, family=my_font),
    strip.background = element_blank(),
    # legend.position = 'none',
    # panel.background = element_blank(),
    # panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )
  p <- p + theme_bw() +  thesis_theme
  lg <- cowplot::get_legend(p + guides(fill = guide_legend(nrow=3, title.position = "left", 
                                                           override.aes = list(shape = 0, size=0.2)))) #nrow=1
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  if(!plotlegend){
    lg_pos <- "none"
    p <- p + ggplot2::theme(legend.position = lg_pos)  
  }else{
    # lg_pos <- "right"
    lg_pos <- "top"
    p <- p + ggplot2::theme(legend.position = lg_pos)  
  }
  
  results <- list(plg=plg, p=p, cols_use=cols_use, df=cls_df, cell_clones=cell_clones)
  # saveRDS(results, paste0(output_dir, datatag, '_',plottitle,"_dlp_prevalence.rds"))
  # nrow <- length(unique(cls_df$timepoint))
  # png(paste0(output_dir,datatag, "_",plottitle, ".png"), height = 2*250, width=2*(50*nrow+20),res = 2*72)
  # print(p)
  # dev.off()
  # 
  # png(paste0(output_dir,datatag, "_",plottitle, "_legend.png"), height = 2*50, 
  #     width=2*(30*length(unique(cls_df$clone_id))+100),res = 2*72)
  # print(plg)
  # dev.off()
  
  return(results)
}
plot_fill_barplot_wholedataset <- function(cell_clones, cols_use, 
                                  output_dir, datatag='SA', plottitle=NULL, plotlegend=F){
  
  # fa='clone_id', xa='timepoint'
  cols_use <- cols_use[unique(cell_clones$clone_id)]
  cls_df <- cell_clones %>%
    dplyr::group_by(clone_id, timepoint, treatment_desc) %>%
    dplyr::summarise(Freq=n())%>%
    as.data.frame()
  # cls_df <- cls_df[gtools::mixedsort(cls_df$clone_id),]
  # dim(cls_df)
  
  cls_df$timepoint <- factor(cls_df$timepoint, levels=gtools::mixedsort(unique(cls_df$timepoint)))
  cls_df$clone_id <- factor(cls_df$clone_id, levels=sort(unique(cls_df$clone_id)))
  ts <- unique(cls_df$treatment_desc)
  ts1 <- ts[grepl('UnRx',ts)]
  ts2 <- ts[grepl('RxH',ts)]
  ts3 <- ts[!ts %in% c(ts1, ts2)]
  cls_df$treatment_desc <- factor(cls_df$treatment_desc, levels=c(ts1, ts3, ts2))
  my_font <- "Helvetica"
  p <- ggplot(cls_df, aes(fill=clone_id, y=Freq, x=timepoint)) + 
    geom_bar(position="fill", stat="identity",width=0.4) + 
    facet_grid(treatment_desc ~ .) + #, space='free', drop=F, scales="free".   . ~ treatment_desc
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) + 
    scale_fill_manual(values = cols_use)#name='Clone '
  p <- p + labs(title=plottitle, x=NULL, y=NULL) #x=xlabel, y=ylabel, y='Clonal fraction'
  
  thesis_theme <- ggplot2::theme(
    # text = element_text(color="black",size = 8, hjust = 0.5, family=my_font),
    # axis.title.x = element_text(color="black",size=8, hjust = 0.5, family=my_font),
    # axis.title.y = element_text(color="black",size=8, hjust = 0.5, family=my_font),
    axis.text.x = element_text(color="black",size=10, hjust = 0.5, family=my_font), #, angle = 90
    axis.text.y = element_text(color="black",size=8, hjust = 0.5, family=my_font), #, angle = 90
    # axis.text.x = element_blank(),
    # axis.ticks.y = element_blank(),
    strip.placement = "outside",
    # axis.line = element_line(colour = "black"),
    # axis.line.x = element_line(colour = "black"),
    # axis.line.y = element_line(colour = "black"),
    # axis.text.y = element_blank(),
    # axis.title = element_blank(),
    plot.title = element_text(color="black",size=11, face="bold",family=my_font, hjust = 0.5),
    legend.title=element_blank(),
    # legend.title=element_text(color="black",size=9, hjust = 0.5, family=my_font),
    legend.text=element_text(color="black",size=9, hjust = 0.5, family=my_font),
    legend.key.height= unit(0.2, 'cm'),
    legend.key.width= unit(0.5, 'cm'),
    strip.text.x = element_text(color="black",size=10, family=my_font),
    strip.text.y = element_text(color="black",size=10, family=my_font),
    strip.background = element_blank(),
    # legend.position = 'none',
    # panel.background = element_blank(),
    # panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )
  p <- p + theme_bw() +  thesis_theme
  lg <- cowplot::get_legend(p + guides(fill = guide_legend(nrow=3, title.position = "left", 
                                                           override.aes = list(shape = 0, size=0.2)))) #nrow=1
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  if(!plotlegend){
    lg_pos <- "none"
    p <- p + ggplot2::theme(legend.position = lg_pos)  
  }else{
    # lg_pos <- "right"
    lg_pos <- "top"
    p <- p + ggplot2::theme(legend.position = lg_pos)  
  }
  
  results <- list(plg=plg, p=p, cols_use=cols_use, df=cls_df, cell_clones=cell_clones)
  # saveRDS(results, paste0(output_dir, datatag, '_',plottitle,"_dlp_prevalence.rds"))
  # nrow <- length(unique(cls_df$timepoint))
  # png(paste0(output_dir,datatag, "_",plottitle, ".png"), height = 2*250, width=2*(50*nrow+20),res = 2*72)
  # print(p)
  # dev.off()
  # 
  # png(paste0(output_dir,datatag, "_",plottitle, "_legend.png"), height = 2*50, 
  #     width=2*(30*length(unique(cls_df$clone_id))+100),res = 2*72)
  # print(plg)
  # dev.off()
  
  return(results)
}  

plot_fill_barplot <- function(cell_clones, cols_use, 
                              plottitle, output_dir, datatag='SA', 
                              plotlegend=F){
  
  # fa='clone_id', xa='timepoint'
  cls_df <- cell_clones %>%
    dplyr::group_by(clone_id, timepoint) %>%
    dplyr::summarise(Freq=n())%>%
    as.data.frame()
  # cls_df <- cls_df[gtools::mixedsort(cls_df$clone_id),]
  # dim(cls_df)
  cls_df$timepoint <- factor(cls_df$timepoint, levels=gtools::mixedsort(unique(cls_df$timepoint)))
  cls_df$clone_id <- factor(cls_df$clone_id, levels=sort(unique(cls_df$clone_id)))
  my_font <- "Helvetica"
  p <- ggplot(cls_df, aes(fill=clone_id, y=Freq, x=timepoint)) + 
    geom_bar(position="fill", stat="identity",width=0.4) + 
    scale_fill_manual(values = cols_use, name='Clone Id: ')
  p <- p + labs(title=plottitle) #x=xlabel, y=ylabel, 
  
  thesis_theme <- ggplot2::theme(
    text = element_text(color="black",size = 8, hjust = 0.5, family=my_font),
    # axis.title.x = element_text(color="black",size=8, hjust = 0.5, family=my_font),
    # axis.title.y = element_text(color="black",size=8, hjust = 0.5, family=my_font),
    axis.text.x = element_text(color="black",size=9, hjust = 0.5, family=my_font), #, angle = 90
    # axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    strip.placement = "outside",
    # axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y= element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(color="black",size=10, face="bold",family=my_font, hjust = 0.5),
    # legend.title=element_text(color="black",size=7, hjust = 0.5, family=my_font),
    legend.text=element_text(color="black",size=9, hjust = 0.5, family=my_font),
    legend.key.height= unit(0.2, 'cm'),
    legend.key.width= unit(0.5, 'cm'),
    strip.text.x = element_text(color="black",size=9, family=my_font),
    strip.text.y = element_text(color="black",size=9, family=my_font),
    # legend.position = 'none',
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )
  p <- p + thesis_theme
  lg <- cowplot::get_legend(p + guides(fill = guide_legend(nrow = 1, title.position = "left", override.aes = list(shape = 0, size=0.6))))
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
  if(!plotlegend){
    lg_pos <- "none"
    p <- p + ggplot2::theme(legend.position = lg_pos)  
  }else{
    lg_pos <- "right"
    p <- p + ggplot2::theme(legend.position = lg_pos)  
  }
  
  results <- list(plg=plg, p=p, cols_use=cols_use, df=cls_df, cell_clones=cell_clones)
  saveRDS(results, paste0(output_dir, datatag, '_',plottitle,"_dlp_prevalence.rds"))
  nrow <- length(unique(cls_df$timepoint))
  png(paste0(output_dir,datatag, "_",plottitle, ".png"), height = 2*250, width=2*(50*nrow+20),res = 2*72)
  print(p)
  dev.off()
  
  png(paste0(output_dir,datatag, "_",plottitle, "_legend.png"), height = 2*50, 
      width=2*(30*length(unique(cls_df$clone_id))+100),res = 2*72)
  print(plg)
  dev.off()
  
  return(results)
}  

plot_all_clones_legends <- function(){
  input_dir <- '/home/htran/storage/datasets/drug_resistance/rna_results/manuscript/umap_figs/figs/'
  
  colorcode <- data.table::fread(paste0(input_dir,'colorcode_total.csv')) %>% as.data.frame()
  
  unassign_clones <- c('unassigned','Unassigned','Un','un','None')
  
  for(pdx in unique(colorcode$datatag)){
    df <- colorcode %>%
      dplyr::filter(datatag==pdx)
    obs_clones <- unique(df$clone_id)
    obs_clones <- obs_clones[!obs_clones %in% unassign_clones]
    df <- df %>%
      dplyr::filter(clone_id %in% obs_clones)
    plot_legends(df, input_dir, pdx)
  }
}

plot_legends <- function(df, output_dir, datatag='SA', legendtitle='Clone Id '){
  my_font <- "Helvetica"
  cols_use <- df$colour
  names(cols_use) <- df$clone_id
  clone_use <- sort(df$clone_id)
  cols_use <- cols_use[clone_use]
  
  df$clone_id <- factor(df$clone_id, levels = clone_use)
  rownames(df) <- df$clone_id
  df <- df[as.character(sort(df$clone_id)),]
  p <- ggplot(df, aes(x=clone_id, y=datatag, color=clone_id)) +
       geom_point(size=8, shape=15) +  #,color=colorcode
       scale_color_manual(values = cols_use,  name=legendtitle) +
       ggplot2::theme(
         legend.title=element_text(color="black",size=9, hjust = 0.5, family=my_font),
         legend.text=element_text(color="black",size=10, hjust = 0.5, family=my_font))
  
  p <- cowplot::get_legend(p + guides(color = guide_legend(ncol=2, override.aes = list(size=4)))) #title.position = "left", 
  plg <- cowplot::ggdraw() + cowplot::draw_plot(p)
  png(paste0(output_dir,datatag, "_dlp_legend.png"), height = 2*30*dim(df)[1], 
      width=2*100,res = 2*72)
  print(plg)
  dev.off()
  saveRDS(list(p=p, plg=plg), paste0(output_dir, datatag,'_clone_legend_plt.rds'))
}

get_library_id <- function(cell_ids) {
  labels <- sapply(strsplit(cell_ids, "-"), function(x) {
    return(x[2])
  })
  return(as.character(labels))
}
get_sample_id <- function(cell_ids) {
  labels <- sapply(strsplit(cell_ids, "-"), function(x) {
    return(x[1])
  })
  return(as.character(labels))
}

get_meta_data <- function(cell_clones_fn, library_grouping_fn, datatag=''){
  # meta_ptx <- data.frame(datatag=c("SA501","SA530","SA604","SA609","SA535","SA1035"), 
  #                        pt=paste0("Pt",rep(1:6,1)))
  meta_pts <- c("SA609"="Pt4","SA535"="Pt5","SA1035"="Pt6",
                "SA501"="Pt1","SA530"="Pt2","SA604"="Pt3")
  cell_clones <- data.table::fread(cell_clones_fn) %>% as.data.frame()
  cell_clones <- cell_clones %>% 
    dplyr::filter(!clone_id %in% c('Un','None','un','unassigned'))
  # unique(cell_clones$clone_id)
  # dim(cell_clones)
  # if(is.null(library_grouping_fn)){
  #   library_grouping_fn <- paste0(results_dir,'library_groupings.csv')
  # }
  grouping_df <- data.table::fread(library_grouping_fn) %>% as.data.frame()
  # colnames(grouping_df)[which(names(grouping_df) == "passage")] <- "timepoint"
  # colnames(grouping_df)[which(names(grouping_df) == "mainsite")] <- "mainsite"
  # colnames(grouping_df)[which(colnames(grouping_df) == "grouping")] <- "library_id"
  # print(colnames(grouping_df))
  cns <- colnames(grouping_df)
  lid <- cns[grepl('library',cns)]
  sid <- cns[grepl('sample',cns)]
  tid <- cns[grepl('treatment',cns)]
  if(datatag=='SA535'){
    lid <- 'grouping'
    grouping_df <- grouping_df %>%
      dplyr::rename(timepoint=passage)
    
    # grouping_df$timepoint <- ifelse(grepl('^X0',grouping_df$timepoint),gsub('0','',grouping_df$timepoint),
    #                                 ifelse(grepl('^UX',grouping_df$timepoint),gsub('U','',grouping_df$timepoint), 
    #                                        grouping_df$timepoint))
    
    grouping_df$timepoint <- paste0('X',gsub('UX|X0|X','',grouping_df$timepoint))
  }
  if(datatag=='SA1035'){
    grouping_df <- grouping_df %>%
      dplyr::rename(timepoint=passage)
  }
  print(unique(grouping_df$timepoint))
  # colnames(grouping_df)[which(colnames(grouping_df) == lid)] <- "library_id"
  # colnames(grouping_df)[which(colnames(grouping_df) == sid)] <- "sample_id"
  # colnames(grouping_df)[which(colnames(grouping_df) == tid)] <- "treatmentSt"
  grouping_df <- grouping_df %>%
    dplyr::rename(library_id=!!sym(lid), sample_id=!!sym(sid), treatmentSt=!!sym(tid))
  # dim(grouping_df)
  # length(unique(grouping_df$treatmentSt))
  unique(grouping_df$sample_id)
  unique(grouping_df$library_id)
  print(colnames(grouping_df))
  cell_clones$library_id <- get_library_id(cell_clones$cell_id)
  cell_clones$sample_id <- get_sample_id(cell_clones$cell_id)
  
  print(dim(cell_clones))
  # unique(cell_clones$treatmentSt)
  cell_clones <- cell_clones %>% left_join(grouping_df, by = c("library_id","sample_id"))
  if(datatag=='SA609'){
    obs_samples <- c('SA609X3XB01584','SA609X4XB03080','SA609X5XB03223',
                     'SA609X6XB03447','SA609X7XB03554','SA609X4XB003083',
                     'SA609X5XB03230','SA609X6XB03404','SA609X7XB03505',
                     'SA609X5XB03231','SA609X6XB03401','SA609X7XB03510')
    # length(obs_samples)
    # sids <- unique(cell_clones$sample_id)
    # unique(grouping_df$sample_id)
    # sum(obs_sids %in% unique(cell_clones$main_sid))
    # sum(obs_samples %in% unique(cell_clones$sample_id))
    # obs_sids <- sapply(obs_samples, function(x){
    #   return(stringr::str_sub(x, nchar(x)-3,nchar(x)))
    # })
    # cell_clones$main_sid <- sapply(cell_clones$sample_id, function(x){
    #   return(stringr::str_sub(x, nchar(x)-3,nchar(x)))
    # })
    # unique(cell_clones$sample)
    cell_clones <- cell_clones %>%
      dplyr::filter(sample_id %in% obs_samples)
    
    # cell_clones2 <- cell_clones %>%
    #   dplyr::filter(clone_id=='E' & treatmentstr!='M')
    # dim(cell_clones2)
    # cell_clones <- dplyr::bind_rows(cell_clones1, cell_clones2)
    print('Replicates excluded!!!')
    # unique(cell_clones$treatmentstr)
    cell_clones <- cell_clones %>% 
      dplyr::filter(clone_id!='A') %>% 
      dplyr::mutate(clone_id=replace(clone_id, clone_id=='R', 'A'))
    treatment_desc <- paste0(meta_pts[datatag],' ', c("Rx","RxH","UnRx"))
    names(treatment_desc) <- c("R","H","U")
    cell_clones$treatment_desc <- treatment_desc[cell_clones$treatmentstr]
    # unique(cell_clones$treatment_desc)
  }else{
    # unique(cell_clones$treatmentSt)
    cell_clones <- cell_clones %>% 
      dplyr::mutate(treatment_desc = case_when(grepl('TU$',treatmentSt) ~ 'RxH',
                                               grepl('T$',treatmentSt) ~ 'Rx',
                                               TRUE  ~ 'UnRx'))
    # cell_clones$treatment_desc <- paste0(meta_pts[datatag],'-', cell_clones$treatment_desc)
    # table(cell_clones$treatmentSt, cell_clones$treatment_desc)
  }
  if(datatag=='SA501'){
    cell_clones$timepoint <- 'X2'
  }
  if(datatag=='SA530'){
    cell_clones$timepoint <- 'X3'
  }
  print(unique(cell_clones$treatment_desc))
  # print(unique(cell_clones$timepoint))
  print(summary(as.factor(cell_clones$clone_id)))
  print(dim(cell_clones))
  
  return(cell_clones)
}
