
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
                        n_genes_to_annotate = 25,
                        n_genes_per_chrom=2,
                        lim_logFC=3.61,
                        min_logFC=0,
                        lim_logFDR=300,
                        sig_FDR=0.01,
                        chroms="all",
                        genetype="all",   # cis, trans or all
                        savecsv=NULL,
                        rel_heights = c(2.5,1)) {
  
  clone_str <- paste(clones, collapse = " vs ")
  
  if (chroms=="all") {
    chrs <- c(as.character(1:22), "X")
  } else {
    chrs=chroms
  }
  
  # include only significant genes
  tt <- tt[tt$FDR <= sig_FDR,]
  

  df_track <- annotables::grch38 %>% 
    dplyr::select(gene_symbol = symbol, chr, start, end) %>% 
    inner_join(tt) %>% 
    #dplyr::filter(chr %in% c(as.character(1:23), "X", "Y")) %>% 
    dplyr::filter(chr %in% chrs) %>% 
    dplyr::mutate(chr = factor(chr, levels = chrs),
                  position = (start + end) / 2)
  
  # MA: sometimes a genes appears more than once, keep only the unique ones
  df_track <- distinct(df_track)
  df_track[df_track$FDR < 10^-{lim_logFDR},]$FDR <- 10^{-runif(length(df_track[df_track$FDR < 10^-{lim_logFDR},]$FDR),lim_logFDR-20,lim_logFDR)}
  
  # Making sure the logFC values are within the limits, otherwise they show grey
  if (length(df_track[df_track$logFC < -lim_logFC,]$logFC) > 0) {
    df_track[df_track$logFC < -lim_logFC,]$logFC <- -lim_logFC
  }
  if (length(df_track[df_track$logFC > lim_logFC,]$logFC) > 0) {
    df_track[df_track$logFC > lim_logFC,]$logFC <- lim_logFC
  }  
  
  # eliminate the genes with logFC < min_logFC
  df_track <- df_track[abs(df_track$logFC) >= min_logFC,]
  
  #df_track$mlog10FDR <- -log10(df_track$FDR)
  #df_track[df_track$mlog10FDR > 290,]$mlog10FDR <- df_track[df_track$mlog10FDR > 290,]$mlog10FDR - runif(1, 10, 20)
  cols <- rev(brewer.pal(n = 7, name = "RdBu") )   # was n=5
  
  # MA: 23 Jan 2021: moving the track plot and gene labels later
  

  # 2 Feb 2021: use grch38! grch37 is missing some genes, such as MARCKS
  #annots <- annotables::grch37
  annots <- annotables::grch38
  
  annots <- dplyr::select(annots, ensembl_gene_id = ensgene, chr, start, end) 
  ## MA: 24 Jan 2021: remove duplicated rows, these would result in spaces in the CN track
  annots <- annots[!duplicated(annots),]
  
  save_df_cnv <- df_cnv
  
  ### MA: include only the genes that have CNV values for both clones
  
  myclones <- unique(df_cnv$cluster)
  genes1 <- df_cnv[df_cnv$cluster==myclones[1],]$ensembl_gene_id
  genes2 <- df_cnv[df_cnv$cluster==myclones[2],]$ensembl_gene_id
  common <- intersect(genes1,genes2)
  df_cnv <- df_cnv[df_cnv$ensembl_gene_id %in% common,]
  
  df_cnv <- inner_join(df_cnv, annots)
  
  df_cnv <- dplyr::select(df_cnv, cluster, ensembl_gene_id, median_cnmode, chr, start, end) %>% 
    group_by(chr) %>% 
    dplyr::mutate(start_order = rank(start)) %>% 
    ungroup() ##%>% 
    ##inner_join(df_cnv)   # not neccesary
  ## MA: because there are 2 conditions that are being compared, start_order is always having a tie
  
  # MA: 24 Jan 2021: rank(start) may give non-integer values because it handles ties using means.
  # Added ties.method="min" -- this doesn't work though
  
  df_cnv$position <- (df_cnv$start+df_cnv$end)/2
  
  if (chroms=="all") {
    chr_levels <- c(as.character(1:22), "X")  
  } else {
    chr_levels <- chroms
  }
  
  df_cnv$chr <- factor(df_cnv$chr, levels = chr_levels)
  
  # removing all the rows that have missing values, for example which are not in our chromosomes of interest
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
  
  ## MA: 16 Jan 2021: setting a factor so that the CN values are plotted in the order given
  df_cnv$cluster <- factor(df_cnv$cluster, levels=c(clones[2],clones[1]))
  
  # If I use position, many genes crowd in one small area
  df_cnv <- dplyr::filter(df_cnv, cluster %in% clones, chr != "Y")
  print(df_cnv)
  
  
  
  
  
  cnv_plot <- 
    ggplot(df_cnv, aes(x = start_order, y = cluster, fill = cnv)) +
    #ggplot(df_cnv, aes(x = position, y = cluster, fill = cnv)) +
    geom_raster() +   ## this works well with start_order
    #geom_tile() +   # this works with position, but all genes crowd in small areas with large spaces in between
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

  ###############

  
  ## MA: creating blank_data for the limits
  blank_data <- NULL
  for (chr in chrs) {
    blank_data <- rbind(blank_data, data.frame(chr=chr, start_order=min(df_cnv[df_cnv$chr==chr,]$start_order), logFC=0))
    blank_data <- rbind(blank_data, data.frame(chr=chr, start_order=max(df_cnv[df_cnv$chr==chr,]$start_order), logFC=0))
  }
  print("Blank data")
  print(blank_data)
  
  df_track$start_order <- NA
  for (geneid in df_track$ensembl_gene_id) {
    df_track[df_track$ensembl_gene_id==geneid,]$start_order <- df_cnv[df_cnv$ensembl_gene_id==geneid,]$start_order[1]  
  }
  df_track <- drop_na(df_track)
  
  # I don't think I need this
  df_track$ciscolor <- "grey"
  df_track[abs(df_track$cnchange) >= 1,]$ciscolor <- "blue"
  
  
  total <- length(df_track$cnchange)
  if (!is.null(savecsv)) {
    write.table(df_track,file=savecsv,sep=",", row.names=FALSE, quote=FALSE)
  }
  ## keep only the cis genes
  if (genetype=="cis") {
    df_track <- df_track[abs(df_track$cnchange) >= 1,]
  } else if (genetype=="trans") {
    df_track <- df_track[abs(df_track$cnchange) < 1,]
  }
  percent <- round(length(df_track$cnchange)/total*100,digits=0)

  sample_name <- paste0(sample_name, " ", percent, "% in-", genetype, " DEGs")
  
  ### Using cis genes and logFC on x
  track_plot <- ggplot(df_track, aes(x = start_order, y = logFC)) +
    geom_point(aes(colour = logFC), size = 1) +   # MA: size was 1
    geom_blank(data = blank_data, aes(x = start_order, y = logFC)) +
    facet_wrap(~ chr, scales = "free_x", 
               strip.position = "bottom",
               # switch = "x", 
               nrow = 1) +
    scale_colour_gradientn(name = glue("logFC {clone_str}"),
                           colours = cols, 
                           values = rescale(c(-2, -0.3, 0, 0.3, 2)),
                           limits = c(-lim_logFC, lim_logFC)) +    
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(limits=c(-lim_logFC,lim_logFC)) +
    theme(strip.background = element_rect(fill = 'white', colour = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=18),  ## was 14      
          strip.placement = "outside",
          legend.position = "top",
          #legend.position = c(0.4,0.8),
          panel.spacing = unit(0.2, 'cm')) +   ## MA: was 0.2
    theme(text = element_text(size = 18)) +    # this was commented out
    geom_hline(yintercept=0, linetype='dashed', col = 'red') +
    labs(x = "Chromosome", y="log 2 fold change", subtitle=sample_name)
    #ggtitle(sample_name)
  
  ## This was before in-cis
  if (FALSE) {
  track_plot <- 
    ggplot(df_track, aes(x = position, y = -log10(FDR))) +
    geom_point(aes(colour = logFC), size = 2) +   # MA: size was 1
    facet_wrap(~ chr, scales = "free_x", 
               strip.position = "bottom",
               # switch = "x", 
               nrow = 1) +
    scale_colour_gradientn(name = glue("log2 fold change, {clone_str}"),
                           colours = cols, 
                           values = rescale(c(-2, -0.3, 0, 0.3, 2)),
                           limits = c(-lim_logFC, lim_logFC)) +
    theme(strip.background = element_rect(fill = 'white', colour = 'white'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=18),  ## was 14      
          strip.placement = "outside",
          legend.position = "top",
          #legend.position = c(0.1,0.7),
          panel.spacing = unit(0.1, 'cm')) +   ## MA: was 0.2
    theme(text = element_text(size = 18)) +    # this was commented out
    ##geom_hline(yintercept=2, linetype='dashed', col = 'red') +
    labs(x = "Chromosome",
         subtitle = sample_name)
    #scale_x_continuous(expand = c(0,0))
  
  }
  
  # MA: This selects the top (by FDR) n_genes_to_annotate genes in the whole genome
  #df_annot <- top_n(df_track, n_genes_to_annotate, -log10(FDR)) 
  # MA: Trying to add a number of genes in each chromosome
  # also make sure we only select significant DEgenes (FDR < 0.01)
  df_annot <- NULL
  #for (chr in chr_levels) {
  #  chr_track <- df_track[df_track$chr==chr & df_track$FDR<0.01 & abs(df_track$cnchange)>=1,]
  #  df_annot <- rbind(df_annot, top_n(chr_track, n_genes_per_chrom, logFC)) 
    
    ## SHowing 2 genes per chromosome
    #chr_track <- df_track[df_track$chr==chr & df_track$FDR<0.01,]
    #df_annot <- rbind(df_annot, top_n(chr_track, n_genes_per_chrom, -log10(FDR)))     
  #}
  
  
  ## label only the genes that have CN change
  
  
  if(!is.null(additional_genes)) {
    df_annot <- df_annot  %>% 
      bind_rows(dplyr::filter(df_track, gene_symbol %in% additional_genes))
  }
  
  ### very good geom_text_repel tutorial
  ## https://ggrepel.slowkow.com/articles/examples.html
  ## thhis is where the genes are added
  track_plot2 <- track_plot +
    ## this makes sure the labels show and do not clip
    coord_cartesian(clip = "off") +
    geom_text_repel(data = df_annot, aes(label = gene_symbol), size = 5,   # was 3.3
                    min.segment.length = 0)
  # point.size=3  default 1, how far from the points to place the label
  ## min.segment.length = 0 always draw the line segments
  
  #not including labels
  main_plot <- cowplot::plot_grid(
    track_plot2 + theme(axis.title.x = element_blank(),
                        strip.text.x = element_blank()),
    cnv_plot,
    ncol = 1,
    rel_heights = rel_heights,
    align = 'v'
  )
  
  main_plot
  #track_plot2
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