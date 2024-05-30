
library(glue)
library(Seurat)
library(RColorBrewer)

# MA: some of the columns can’t be labeled certain keywords
# This function relabels some of the columns in a SingleCellExperiment object so it can be properly subset and used
fixsce <- function (sce) {
  # some of them don;t have these, removing them for now
  # rowData(sce)$chromosome_name <- rowData(sce)$chr
  # rowData(sce)$start_position <- rowData(sce)$start
  # rowData(sce)$end_position <- rowData(sce)$end
  drops <- c("chr","start","end","strand","symbol","ensgene","entrez","biotype","description")
  rdcopy <- rowData(sce)
  rowData(sce) <- rdcopy[ , !(names(rdcopy) %in% drops)]
  return(sce)
}

#' This function makes a plot of the log fold change between to clones by chromosomal
#' location, as well as the copy number profiles of the two clones
#' 
#' @param tt Differential expression results (edgeR)
#' @param df_cnv CNV data as dataframe
#' @param clones Clones under consideration (length 2)
#' @param sample_name Name of sample / passage
plot_gex_cnv <- function(tt,
                        df_cnv,
                        clones,
                        sample_name,
                        additional_genes = NULL,
                        n_genes_to_annotate = 25) {
  
  clone_str <- paste(clones, collapse = " vs ")
  
  chrs <- c(as.character(1:23), "X", "Y")
  df_track <- annotables::grch37 %>% 
    dplyr::select(gene_symbol = symbol, chr, start, end) %>% 
    inner_join(tt) %>% 
    dplyr::filter(chr %in% c(as.character(1:23), "X", "Y")) %>% 
    dplyr::mutate(chr = factor(chr, levels = chrs),
                  position = (start + end) / 2)
  
  # MA: sometimes a genes appears more than once, keep only the unique ones
  df_track <- distinct(df_track)
  
  cols <- rev(brewer.pal(n = 5, name = "RdBu") )
  
  track_plot <- ggplot(df_track, aes(x = position, y = -log10(FDR))) +
    geom_point(aes(colour = logFC), size = 2) +   # MA: size was 1
    facet_wrap(~ chr, scales = "free_x", 
               strip.position = "bottom",
               # switch = "x", 
               nrow = 1) +
    scale_colour_gradientn(name = glue("log2 fold change, {clone_str}"),
                           colours = cols, 
                           values = rescale(c(-2, -0.3, 0, 0.3, 2)),
                           limits = c(-3.61, 3.61)) +
    theme(strip.background = element_rect(fill = 'white', colour = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=14),        
          strip.placement = "outside",
          legend.position = "top",
          panel.spacing = unit(0.1, 'cm')) +   ## MA: was 0.2
    theme(text = element_text(size = 18)) +
    labs(x = "Chromosome",
         subtitle = sample_name) +
    scale_x_continuous(expand = c(0,0))
  
  df_annot <- top_n(df_track, n_genes_to_annotate, -log10(FDR)) 
  
  
  print("The labeled genes")
  print(df_annot)  
  
  if(!is.null(additional_genes)) {
    df_annot <- df_annot  %>% 
      bind_rows(dplyr::filter(df_track, gene_symbol %in% additional_genes))
  }
  
  ## thhis is where the genes are added
  track_plot2 <- track_plot +
    geom_text_repel(data = df_annot, aes(label = gene_symbol), size = 3.3)
  
  annots <- annotables::grch37
  
  annots <- dplyr::select(annots, ensembl_gene_id = ensgene, chr, start, end) 
  
  df_cnv <- inner_join(df_cnv, annots)
  
  df_cnv <- dplyr::select(df_cnv, ensembl_gene_id, chr, start) %>% 
    group_by(chr) %>% 
    dplyr::mutate(start_order = rank(start)) %>% 
    ungroup() %>% 
    inner_join(df_cnv)
  
  chr_levels <- c(as.character(1:23), "X", "Y")
  
  df_cnv$chr <- factor(df_cnv$chr, levels = chr_levels)
  
  df_cnv <- drop_na(df_cnv)
  
  # MA: 11 Apr 2020, setting the copy number colors that were used in the heatmap
  # TO COME BACK
  #print("data frame")
  #cnv_colors <- c('#4880B8', '#A7C9DF','#CCCCCC','#F5CE93','#ED9364','#D2553E','#A42116','#8B1A43','#CB3576','#D06CAD','#C196C4','#D0BAD8')
  #cnv_cols <- data.frame(cn=0:11,color<-cnv_cols)
  #colnames(cnv_cols) <- c("cn","color")
  #levels(cnv_cols$color) <- cnv_colors
  
  cnv_colors <- c('#4880B8', '#A7C9DF','#CCCCCC','#F5CE93','#ED9364','#D2553E','#A42116','#8B1A43','#CB3576','#D06CAD','#C196C4','#D0BAD8')

  cnv_cols <- c('0'='#4880B8', '1'='#A7C9DF','2'='#CCCCCC','3'='#F5CE93','4'='#ED9364',
                '5'='#D2553E','6'='#A42116','7'='#8B1A43','8'='#CB3576','9'='#D06CAD','10'='#C196C4','11'='#D0BAD8')
  #levels(cnv_cols) <- 0:11
#  cnv_cols <- c("0" = "#2166ac",
#                "1" = "#92c5de", 
#                "2" = "grey80", 
#                "3" = "#f4a582", 
#                "4" = "#d6604d",
#                "5" = "#b2182b",
#                "6+" = "#67001f")
  
  df_cnv$cnv <- as.character(round(df_cnv$median_cnmode))
  levels(df_cnv$cnv) <- 0:(length(cnv_colors)-1)
  #print(levels(df_cnv$cnv))
  
  # MA: removing the 6+ restriction (clonealign still has that restriction though)
  #df_cnv$cnv[round(df_cnv$median_cnmode) >= 6] <- "6+"
  
  cnv_plot <- dplyr::filter(df_cnv, cluster %in% clones, chr != "Y") %>% 
    ggplot(aes(x = start_order, y = cluster, fill = cnv)) +
    geom_raster() +
    facet_wrap(~ chr, scales = "free_x", nrow = 1, switch = "x") +
    #theme(legend.position = "top", axis.text.x = element_blank()) +
    scale_y_discrete(expand = c(0, 0)) +
    #scale_fill_manual(values=cnv_colors, name = "Copy number", guide = 'legend',labels = 0:(length(cnv_colors)-1),drop=FALSE)  +
    #          theme(legend.position = "bottom") +   
    scale_fill_manual(values = cnv_cols, name = "Copy number") +  #, labels = 0:(length(cnv_cols)-1),drop=FALSE) +
    labs(x = "chromosome") +
    theme(strip.background = element_rect(fill = 'white', colour = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=14), 
          strip.placement = "outside",
          legend.position = "none",
          panel.spacing = unit(c(0.1), 'cm')) +   # MA: was 0.2
    theme(text = element_text(size = 18)) +
    labs(x = "Chromosome", y = "Clone")+
    guides(fill = guide_legend(nrow = 1)) +
    scale_x_continuous(expand = c(0,0))

  
  main_plot <- cowplot::plot_grid(
    track_plot2 + theme(axis.title.x = element_blank(),
                        strip.text.x = element_blank()),
    cnv_plot,
    ncol = 1,
    rel_heights = c(2.5,1),
    align = 'v'
  )
  
  main_plot
}

#' Get cluster colours
get_cluster_colours <- function(nClusters=9) {
  if (nClusters > 8) {
    clust_colours <- colorRampPalette(brewer.pal(8, "Set2"))(nClusters)
  } else {
    clust_colours <- brewer.pal(nClusters, "Set2")
  }
  
  if (nClusters == 9) {
    clust_colours[[5]] <- '#564527'
    
    # change I to pink too
    clust_colours[[9]] <- '#FFC0CB'
  }
  
  names(clust_colours) <- LETTERS[seq_along(clust_colours)]
  clust_colours
} 

#' Basic seurat clustering
seurat_cluster <- function(sce, 
                           resolution = 0.8,
                           dims = 1:10,
                           reduction = "mnn",
                           algorithm = 1,
                           return_SCE = TRUE) {
  
  seu <- as.seurat(sce)
  # seu <- FindNeighbors(seu, dims = dims, reduction = reduction)
  seu <- FindClusters(seu, resolution = resolution, algorithm = algorithm)
  
  sce2 <- as.SingleCellExperiment(seu)
  sce$seurat_clusters <- sce2$seurat_clusters
  if(return_SCE) {
    return(sce)
  } else {
    return(sce$seurat_clusters)
  }
}