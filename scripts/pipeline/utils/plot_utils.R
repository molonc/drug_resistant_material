suppressPackageStartupMessages({
  library(tidyverse)
  library(annotables)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  require(scales)
  require(ggrepel)
  require(stringr)
})
library(extrafont)
font_import(prompt=F, paths ='/usr/share/fonts/truetype/myfonts/') # import Helvetica font
fonts()
# devtools::install_github("stephenturner/annotables")
get_top_genes <- function(de_genes, ntop=NULL){
  # de_genes <- de_genes %>%
  #   dplyr::filter(logFC>minLogFC & FDR<FDR_cutoff & PValue<pValueThrs)
  de_genes <- de_genes[order(de_genes$logFC,-log10(de_genes$FDR), decreasing = T),] 
  if(ntop > nrow(de_genes)){
    ntop = nrow(de_genes)
  }
  if(!is.null(ntop)){
    return(de_genes[1:ntop,])  
  }else{  # select all resistant genes with logFC > 0
    return(de_genes)
  }
  
}

plot_gex_cnv_v2 <- function(deg_fn,
                         df_cnv,
                         clones,
                         save_dir,
                         sample_name='',
                         clone_str=NULL,
                         additional_genes_fn = NULL,
                         n_genes_to_annotate = 10, 
                         plttitle=NULL,
                         gene_type='in-cis') {
  print(deg_fn)
 
  if(is.null(plttitle)){
    tl <- element_blank()
  }else{
    tl <- element_text(size=15, colour = "black", face='bold')  # 
  }
  if(is.null(clone_str)){
    clone_str <- paste(clones, collapse = " vs ")
  }
  if(is.null(additional_genes_fn)){
    additional_genes <- NULL
  }else{
    print(additional_genes_fn)
    additional_genes_df <- read.csv(additional_genes_fn, stringsAsFactors = F, check.names = F)
    additional_genes <- unique(additional_genes_df[,1])
    additional_genes <- gsub(' ','', additional_genes)
    print(additional_genes)
  }
  
  deg_df <- read.csv(deg_fn, check.names=F, stringsAsFactors=F) #, row.names = 1
  print(dim(deg_df))
  # length(intersect(deg_df$ensembl_gene_id, pcnv_609$df_cnv$ensembl_gene_id))
  # deg_df <- deg_df %>%
  #   dplyr::filter(abs(logFC)>0.25 & FDR<0.01 & PValue<0.05)
  # deg_df <- deg_df %>%
  #   dplyr::filter(abs(logFC)>0.5 & FDR<0.01 & PValue<0.05)
  # print(dim(deg_df))
  print(gene_type)
  if(gene_type=='in-cis'){
    deg_df <- deg_df %>%
      dplyr::filter(!is.na(classified_gene_dlp))
  }else{
    deg_df <- deg_df %>%
      dplyr::filter(is.na(classified_gene_dlp))
  }
  print(dim(deg_df))
  # deg_df <- deg_df[abs(deg_df$logFC)>0.25 & deg_df$PValue<0.05 & deg_df$FDR<0.01,]
  # deg_df <- deg_df %>% rownames_to_column("ensembl_gene_id")
  # colnames(deg_df)[which(names(deg_df) == "gene_symb")] <- "gene_symbol"
  # summary(as.factor(!is.na(deg_df$classified_gene_dlp)))
  # length(unique(df_cnv$ensembl_gene_id))
  # colnames(df_cnv)
  # df_cnv1 <- df_cnv
  # df_cnv <- df_cnv %>%
  #   dplyr::group_by(ensembl_gene_id) %>%
  #   dplyr::summarise(is_variant=var(cnv)>0) %>%
  #   dplyr::filter(is_variant==TRUE)
  
  # deg_df$gene_type <- ifelse(deg_df$ensembl_gene_id %in% unique(df_cnv$ensembl_gene_id),'in-cis','in-trans')
  # print(summary(as.factor(deg_df$gene_type)))
  # deg_df$is_genome <- !is.na(deg_df$classified_gene_dlp)
  # print(summary(as.factor(deg_df$classified_gene)))
 
  chrs <- c(as.character(1:22), "X")
  # annots <- dplyr::select(annots, ensembl_gene_id = ensgene, chr, start, end) 
  annots <- annotables::grch38
  annots <- dplyr::select(annots, ensembl_gene_id = ensgene, chr, start, end) 
  ## MA: 24 Jan 2021: remove duplicated rows, these would result in spaces in the CN track
  annots <- annots[!duplicated(annots),]
  df_track <- deg_df %>%
    inner_join(annots, by=c('ensembl_gene_id')) %>%
    dplyr::filter(chr %in% chrs) %>%
    dplyr::mutate(position = (start + end) / 2)#chr = factor(chr, levels = chrs),
  # unique(df_track$chr)
  # df_track <- df_track[!is.na(df_track$chr),]
  # df_track <- annotables::grch38 %>%
  #   # dplyr::select(gene_symbol = symbol, chr, start, end) %>%
  #   dplyr::select(ensembl_gene_id = ensgene, chr, start, end) %>%
  #   inner_join(deg_df) %>%
  #   dplyr::filter(chr %in% chrs) %>%
  #   dplyr::mutate(chr = factor(chr, levels = chrs),
  #                 position = (start + end) / 2)
  
  # MA: sometimes a genes appears more than once, keep only the unique ones
  df_track <- distinct(df_track)
  
  # obs_chrs <- c(as.character(13),as.character(22))
  # df_track_obs <- df_track %>%
  #                 dplyr::filter(chr %in% obs_chrs)
  
  
  cols <- rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu"))
  # cols <- c("#021621", "#92C5DE", "#F7F7F7", "#F4A582", "#610614")
  #blue: #053885, red: #8f0319
  # df_track$mlog10FDR <- -log10(df_track$FDR)
  # df_track$mlog10FDR <- sapply(df_track$mlog10FDR, function(x) replace(x, is.infinite(x), 270))
  # df_track$mlog10FDR <- sapply(df_track$mlog10FDR, function(x) replace(x, x>270, 270))
  
  # Rescale data for plotting
  maxLogFC <- 3.5
  df_track$logFC <- sapply(df_track$logFC, function(x) replace(x, x > maxLogFC, maxLogFC))
  df_track$logFC <- sapply(df_track$logFC, function(x) replace(x, x < (-maxLogFC), -maxLogFC))
  # summary(df_track$logFC)
  # df_track <- drop_na(df_track)
  df_track$chr <- as.character(df_track$chr)
  df_track$chr <- factor(df_track$chr, levels = chrs)
  blank_data <- NULL
  if(gene_type=='in-cis'){
    for (chr in chrs) {
      blank_data <- rbind(blank_data, data.frame(chr=chr, start_order=min(df_cnv[df_cnv$chr==chr,'start_order']), logFC=0))
      blank_data <- rbind(blank_data, data.frame(chr=chr, start_order=max(df_cnv[df_cnv$chr==chr,'start_order']), logFC=0))
    }
    df_track$start_order <- 0.0
    genes_use <- intersect(df_track$ensembl_gene_id, df_cnv$ensembl_gene_id)
    # length(genes_use)
    # rownames(df_track) <- df_track$ensembl_gene_id
    # rownames(df_cnv) <- df_cnv$ensembl_gene_id
    # length(unique(df_track$ensembl_gene_id))
    for (geneid in genes_use) {
      t <- df_cnv[df_cnv$ensembl_gene_id==geneid,]$start_order
      df_track[df_track$ensembl_gene_id==geneid,'start_order'] <- as.numeric(t[1])
      # df_track[geneid,'start_order'] <- as.numeric(t[1])
    }
    track_plot <- ggplot(df_track, aes(x = start_order, y = logFC)) +
      geom_point(aes(colour = logFC), size = 1.3) +
      geom_blank(data = blank_data, aes(x = start_order, y = logFC)) 
  } else{
    df_cnv <- df_cnv %>%
      dplyr::mutate(position = (start + end) / 2)
    for (chr in chrs) {
      blank_data <- rbind(blank_data, data.frame(chr=chr, position=min(df_cnv[df_cnv$chr==chr,'position']), logFC=0))
      blank_data <- rbind(blank_data, data.frame(chr=chr, position=max(df_cnv[df_cnv$chr==chr,'position']), logFC=0))
    }
    track_plot <- ggplot(df_track, aes(x = position, y = logFC)) +
      geom_point(aes(colour = logFC), size = 1.3) +
      geom_blank(data = blank_data, aes(x = position, y = logFC)) 
  }
  track_plot <- track_plot + 
    geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5) + 
    facet_wrap(~ factor(chr, levels = chrs), scales = "free_x",
               strip.position = "bottom",
               switch = "x",
               nrow = 1, drop = F) +
    # facet_grid(~ chr, scales = "free_x", space='free',
    #            # strip.position = "bottom",
    #            switch = "x"
    #            # nrow = 1
    # ) +
    scale_colour_gradientn(name = 'log2FC ',  #paste0(clone_str,'  log2FC ')
                           colours = cols, 
                           values = rescale(c(-maxLogFC, -1, 0, 1, maxLogFC)), #rescale(c(-2, -0.3, 0, 0.3, 2)
                           limits = c(-maxLogFC, maxLogFC)) 
    # scale_shape_manual(name= 'gene type ',values=c(17, 16))+
    # theme(strip.background = element_rect(fill = 'white', colour = 'white'),
    #       plot.title = tl,
    #       plot.subtitle = element_text(size=13, colour = "black", face='bold'),  # 
    #       axis.text.x = element_blank(),
    #       axis.ticks.x = element_blank(),
    #       axis.text.y = element_text(size=12, colour = "black", hjust = 0.5),   
    #       axis.title.y = element_text(size=12, colour = "black", hjust = 0.5),   
    #       axis.line = element_line(colour = "black"),
    #       strip.placement = "outside",
    #       legend.position = "top",
    #       legend.title = element_text(color="black", size=13, hjust = 0.5), 
    #       legend.text = element_text(color="black", size=10, hjust = 0.5), 
    #       panel.background = element_blank(), 
    #       panel.spacing = unit(0.1, 'cm')) +   ## MA: was 0.2
    # theme(text = element_text(size = 18)) +
    my_font <- "Helvetica"
    
    # lg_pos <- "none"
    lg_pos <- "top"
    thesis_theme <- ggplot2::theme(
      text = element_text(color="black",size = 8, hjust = 0.5, family=my_font),
      axis.title.x = element_text(color="black",size=8, hjust = 0.5, family=my_font),
      axis.title.y = element_text(color="black",size=8, hjust = 0.5, family=my_font),
      # axis.text.x = element_text(color="black",size=7, hjust = 0.5, family=my_font, angle = 90),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.placement = "outside",
      axis.line = element_line(colour = "black"),
      axis.text.y = element_text(color="black",size=7, hjust = 0.5, family=my_font),
      plot.title = element_text(color="black",size=10, face="bold", hjust=0, family=my_font),
      legend.title=element_text(color="black",size=7, hjust = 0.5, family=my_font),
      legend.text=element_text(color="black",size=7, hjust = 0.5, family=my_font),
      strip.text.x = element_text(color="black",size=9, family=my_font),
      strip.text.y = element_text(color="black",size=9, family=my_font),
      legend.spacing.x = unit(0.1, 'mm'),
      legend.spacing.y = unit(0.1, 'mm'),
      legend.key.height=unit(0.3,"cm"), # line
      # legend.key.height = unit(1, 'cm'), #change legend key height
      legend.key.width = unit(0.6, 'cm'), #change legend key width
      legend.position = lg_pos,
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-5,-5,-5,-5),
      panel.background = element_blank(), 
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    )
   track_plot <- track_plot + thesis_theme
   track_plot <- track_plot + labs(x = "Chromosome", y="log2 FC", title = plttitle) +  #subtitle = sample_name
                              scale_x_continuous(expand = c(0,0)) #+ 
                              #guides(colour = guide_legend(nrow = 1, override.aes = list(size=0.05)))
  
  # track_plot
  # df_annot <- top_n(df_track, n_genes_to_annotate, -log10(FDR))
  # df_annot <- top_n(df_track, n_genes_to_annotate, logFC) 
  # df_track <- df_track[order(df_track$logFC, df_track$mlog10FDR, decreasing = T),] 
  df_track <- df_track[order(abs(df_track$logFC), decreasing = T),] 
  df_annot <- df_track[1:n_genes_to_annotate,]
  print(df_annot$gene_symbol)
  # dim(df_annot)
  
  print("The labeled genes")
  if(!is.null(additional_genes)) {
    df_annot <- df_annot  %>% 
      bind_rows(dplyr::filter(df_track, gene_symbol %in% additional_genes))
  }
  df_annot <- df_annot %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::distinct()
  print(dim(df_annot))
  
  ## thhis is where the genes are added
  # track_plot2 <- track_plot +
  #   geom_text_repel(data = df_annot, aes(label = gene_symbol), 
  #                   size = 2.8)#, xlim = c(NA, Inf),ylim = c(-Inf, Inf),box.padding = 0.2
  # 
  ### very good geom_text_repel tutorial
  ## https://ggrepel.slowkow.com/articles/examples.html
  ## thhis is where the genes are added
  track_plot2 <- track_plot +
    ## this makes sure the labels show and do not clip
    coord_cartesian(clip = "off") +
    geom_text_repel(data = df_annot, aes(label = gene_symbol), size = 3.5,   # was 3.3
                    min.segment.length = 0, max.overlaps = Inf,
                    color='black', segment.alpha = 0.2,
                    family = my_font)
  
  
  return(track_plot2)
}


#' This function makes a plot of the log fold change between to clones by chromosomal
#' location, as well as the copy number profiles of the two clones
#' 
#' @param tt Differential expression results (edgeR)
#' @param df_cnv CNV data as dataframe
#' @param clones Clones under consideration (length 2)
#' @param sample_name Name of sample / passage
plot_gex_cnv <- function(deg_fn,
                         df_cnv,
                         cnv_plot,
                         clones,
                         save_dir,
                         sample_name='',
                         clone_str=NULL,
                         additional_genes_fn = NULL,
                         n_genes_to_annotate = 20) {
  print(deg_fn)
  if(is.null(clone_str)){
    clone_str <- paste(clones, collapse = " vs ")
  }
  if(is.null(additional_genes_fn)){
    additional_genes <- NULL
  }else{
    print(additional_genes_fn)
    additional_genes_df <- read.csv(additional_genes_fn, stringsAsFactors = F, check.names = F)
    additional_genes <- additional_genes_df[,1]
    additional_genes <- gsub(' ','', additional_genes)
    print(additional_genes)
  }
  
  deg_df <- read.csv(deg_fn, check.names=F, stringsAsFactors=F) #, row.names = 1
  print(dim(deg_df))
  deg_df <- deg_df %>%
    dplyr::filter(abs(logFC)>0.25 & FDR<0.01 & PValue<0.05)
  print(dim(deg_df))
  # deg_df <- deg_df[abs(deg_df$logFC)>0.25 & deg_df$PValue<0.05 & deg_df$FDR<0.01,]
  # deg_df <- deg_df %>% rownames_to_column("ensembl_gene_id")
  # colnames(deg_df)[which(names(deg_df) == "gene_symb")] <- "gene_symbol"
  # summary(as.factor(!is.na(deg_df$classified_gene_dlp)))
  # length(unique(df_cnv$ensembl_gene_id))
  # colnames(df_cnv)
  # df_cnv1 <- df_cnv
  df_cnv <- df_cnv %>%
            dplyr::group_by(ensembl_gene_id) %>%
            dplyr::summarise(is_variant=var(cnv)>0) %>%
            dplyr::filter(is_variant==TRUE)
  
  deg_df$gene_type <- ifelse(deg_df$ensembl_gene_id %in% unique(df_cnv$ensembl_gene_id),'in-cis','in-trans')
  print(summary(as.factor(deg_df$gene_type)))
  # deg_df$is_genome <- !is.na(deg_df$classified_gene_dlp)
  # print(summary(as.factor(deg_df$classified_gene)))
  
  chrs <- c(as.character(1:22), "X")
  # annots <- dplyr::select(annots, ensembl_gene_id = ensgene, chr, start, end) 
  df_track <- annotables::grch38 %>% 
    # dplyr::select(gene_symbol = symbol, chr, start, end) %>% 
    dplyr::select(ensembl_gene_id = ensgene, chr, start, end) %>% 
    inner_join(deg_df) %>% 
    dplyr::filter(chr %in% chrs) %>% 
    dplyr::mutate(chr = factor(chr, levels = chrs),
                  position = (start + end) / 2)
  
  # MA: sometimes a genes appears more than once, keep only the unique ones
  df_track <- distinct(df_track)
  
  # obs_chrs <- c(as.character(13),as.character(22))
  # df_track_obs <- df_track %>%
  #                 dplyr::filter(chr %in% obs_chrs)
  
  
  cols <- rev(brewer.pal(n = 5, name = "RdBu"))
  df_track$mlog10FDR <- -log10(df_track$FDR)
  df_track$mlog10FDR <- sapply(df_track$mlog10FDR, function(x) replace(x, is.infinite(x), 270))
  df_track$mlog10FDR <- sapply(df_track$mlog10FDR, function(x) replace(x, x>270, 270))
  
  df_track$logFC <- sapply(df_track$logFC, function(x) replace(x, x > 3.6, 3.6))
  df_track$logFC <- sapply(df_track$logFC, function(x) replace(x, x < (-3.6), -3.6))
  
  
  track_plot <- ggplot(df_track, aes(x = position, y = mlog10FDR)) +  #-log10(FDR)
    geom_point(aes(colour = logFC, shape=gene_type), size = 1.8) +   # MA: size was 1
    facet_wrap(~ chr, scales = "free_x",
               strip.position = "bottom",
               switch = "x",
               nrow = 1) +
    # facet_grid(~ chr, scales = "free_x", space='free',
    #            # strip.position = "bottom",
    #            switch = "x"
    #            # nrow = 1
    # ) +
    scale_colour_gradientn(name = paste0(clone_str,', log2 FC'),
                           colours = cols, 
                           values = rescale(c(-2, -0.3, 0, 0.3, 2)), #rescale(c(-2, -0.3, 0, 0.3, 2)
                           limits = c(-3.6, 3.6)) +
    scale_shape_manual(name= 'gene type ',values=c(17, 16))+
    theme(strip.background = element_rect(fill = 'white', colour = 'white'),
          plot.subtitle = element_text(size=14, colour = "black", face='bold'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=12, colour = "black"),   
          axis.line = element_line(colour = "black"),
          strip.placement = "outside",
          legend.position = "top",
          legend.title = element_text(color="black", size=14, hjust = 0.5), 
          legend.text = element_text(color="black", size=10, hjust = 0.5), 
          panel.background = element_blank(), 
          panel.spacing = unit(0.1, 'cm')) +   ## MA: was 0.2
    theme(text = element_text(size = 18)) +
    labs(x = "Chromosome", y="-log10(FDR)",
         subtitle = sample_name) +
    scale_x_continuous(expand = c(0,0))
  
  # track_plot
  # df_annot <- top_n(df_track, n_genes_to_annotate, -log10(FDR))
  # df_annot <- top_n(df_track, n_genes_to_annotate, logFC) 
  df_track <- df_track[order(df_track$logFC, df_track$mlog10FDR, decreasing = T),] 
  df_annot <- df_track[1:n_genes_to_annotate,]
  print(df_annot$gene_symbol)
  # dim(df_annot)
  
  print("The labeled genes")
  # print(df_annot)  
  
  if(!is.null(additional_genes)) {
    df_annot <- df_annot  %>% 
      bind_rows(dplyr::filter(df_track, gene_symbol %in% additional_genes))
  }
  df_annot <- df_annot %>%
               dplyr::group_by(gene_symbol) %>%
               dplyr::distinct()
  print(dim(df_annot))
  
  ## thhis is where the genes are added
  track_plot2 <- track_plot +
    geom_text_repel(data = df_annot, aes(label = gene_symbol), size = 3.7)
  
  # return(track_plot2)
  main_plot <- cowplot::plot_grid(
    track_plot2 + theme(axis.title.x = element_blank(),
                        strip.text.x = element_blank()),
    cnv_plot,
    ncol = 1,
    rel_heights = c(2.5,1),
    align = 'v'
  )
  save_fn <- paste0(paste(sample_name, collapse = '_'),'_',paste(clone_str, collapse = '_'))
  png(paste0(save_dir,save_fn,"_track_plot.png"), height = 2*500, width=2*1300, res = 2*72)
  print(main_plot)
  dev.off()
  
  # write.csv(df_track_obs, file=paste0(save_dir,save_fn,"_chr_13_22.csv"), quote=F, row.names = F)
  saveRDS(main_plot, paste0(save_dir,save_fn,"_track_plot.rds"))
  
}


plot_DE_genes_chr <- function(save_dir, cnv_mat, input_deg_fn, obs_clones=c('E','H'), 
                              clone_str='',datatag='',
                              additional_genes = NULL, n_genes_to_annotate = 40,
                              outlier_thres=4.25){
  
  # rownames(cnv_mat) <- cnv_mat$ensembl_gene_id
  # cnv_mat$ensembl_gene_id <- rownames(cnv_mat)
  deg_df <- read.csv(input_deg_fn, check.names=F, stringsAsFactors=F) #, row.names = 1
  print(dim(deg_df))
  deg_df <- deg_df[abs(deg_df$logFC)>0.25 & deg_df$PValue<0.05 & deg_df$FDR<0.01,]
  # deg_df <- deg_df[deg_df$logFC<=outlier_thres,]
  # summary(abs(deg_df$p_adj_val))
  # View(head(deg_df))
  # deg_df <- deg_df %>% rownames_to_column("ensembl_gene_id")
  # colnames(deg_df)[which(names(deg_df) == "gene_symb")] <- "gene_symbol"
  summary(as.factor(!is.na(deg_df$classified_gene_dlp)))
  # deg_df$is_genome <- deg_df$ensembl_gene_id %in% cnv_mat$ensembl_gene_id
  deg_df$is_genome <- !is.na(deg_df$classified_gene_dlp)
  print(summary(as.factor(deg_df$classified_gene)))
  # deg_df$classified_gene[1:3]
  # deg_df$variance <- ifelse(deg_df$is_genome==TRUE,'incis',
  #                           ifelse(deg_df$is_genome==FALSE,'intrans','unknown'))
  deg_df$variance <- deg_df$classified_gene
  print(summary(as.factor(deg_df$variance)))
  if(clone_str==''){
    clone_str <- paste(obs_clones, collapse = " vs ")
    clone_str <- paste0(datatag, ':  ', clone_str)
    print(clone_str)
    # clone_str <- 'SA609:  R-UTTT vs H-UUUU'
  }
  
  
  # chrs <- c(as.character(1:22), "X", "Y")
  chrs <- c(as.character(1:22))
  df_track <- annotables::grch37 %>% 
    dplyr::select(gene_symbol = symbol, chr, start, end) %>% 
    inner_join(deg_df) %>% 
    dplyr::filter(chr %in% chrs) %>% 
    dplyr::mutate(chr = factor(chr, levels = chrs),
                  position = (start + end) / 2)
  
  # MA: sometimes a genes appears more than once, keep only the unique ones
  df_track <- distinct(df_track)
  # dim(df_track)
  cols <- rev(brewer.pal(n = 3, name = "RdBu") )
  cols <- cols[1:length(unique(deg_df$is_genome))]
  # cols <- c("#E82503", "#093CB4")
  track_plot <- ggplot(df_track, aes(x = position, y = logFC)) +  #-log10(p_val_adj)
    geom_point(aes(colour = variance), size = 1.3) +   # MA: size was 1
    facet_grid(~ chr, scales = "free_x", space='free',
               # strip.position = "bottom",
               switch = "x"
               # nrow = 1
    ) +
    # scale_color_manual(values = cols) + 
    # scale_colour_gradientn(name = glue("is variant genome"),
    #                        colours = cols, 
    #                        values = c('TRUE','FALSE')) +
    theme(#strip.background = element_blank(),
      strip.background = element_rect(fill = 'white', colour = 'white'),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=10),   
      legend.text = element_text(size=10),
      legend.title = element_text(size=10),
      strip.placement = "outside",
      legend.position = "top",
      panel.spacing = unit(0.1, 'cm'),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # panel.background = element_blank(), 
      axis.line = element_line(colour = "black")) +   ## MA: was 0.2
    theme(text = element_text(size = 18)) +
    # theme_classic() +
    labs(x = "Chromosome",
         subtitle = clone_str) +
    scale_x_continuous(expand = c(0,0))
  
  
  
  df_annot <- dplyr::top_n(df_track[df_track$is_genome,], n_genes_to_annotate, abs(logFC))
  
  
  print("The labeled genes")
  print(dim(df_annot))
  
  if(!is.null(additional_genes)) {
    df_annot <- df_annot  %>% 
      bind_rows(dplyr::filter(df_track, gene_symbol %in% additional_genes))
  }
  
  ## this is where the genes are added
  track_plot2 <- track_plot +
    geom_text_repel(data = df_annot, aes(label = gene_symbol), size = 2)
  
  # base_name <- paste(obs_clones, collapse = "_vs_")
  # png(paste0(save_dir,"DE_chr_",base_name,".png"), height = 2*450, width=2*1200,res = 2*72)
  # print(track_plot2)
  # dev.off()
  
  return(track_plot2)
}


plot_DE_genes_chr_logFC <- function(save_dir, cnv_mat, input_deg_fn, obs_clones=c('E','H'), 
                              sample_name='',datatag='',
                              additional_genes = NULL, n_genes_to_annotate = 60,
                              outlier_thres=4.25){
  
  rownames(cnv_mat) <- cnv_mat$ensembl_gene_id
  deg_df <- read.table(input_deg_fn, check.names=F, stringsAsFactors=F, header=T, row.names = 1)
  print(dim(deg_df))
  deg_df <- deg_df[deg_df$logFC<=outlier_thres,]
  
  deg_df$mlog10FDR <- -log10(deg_df$FDR)
  deg_df$mlog10FDR <- sapply(deg_df$mlog10FDR, function(x) replace(x, is.infinite(x), 300))
  # summary(abs(deg_df$p_adj_val))
  # View(head(deg_df))
  # deg_df <- deg_df %>% rownames_to_column("ensembl_gene_id")
  # colnames(deg_df)[which(names(deg_df) == "gene_symb")] <- "gene_symbol"
  deg_df$is_genome <- deg_df$ensembl_gene_id %in% cnv_mat$ensembl_gene_id
  print(summary(as.factor(deg_df$is_genome)))
  deg_df$variance <- ifelse(deg_df$is_genome==TRUE,'incis',
                            ifelse(deg_df$is_genome==FALSE,'intrans','unknown'))
  
  print(summary(as.factor(deg_df$variance)))
  
  clone_str <- paste(obs_clones, collapse = " vs ")
  clone_str <- paste0(datatag, ':  ', clone_str)
  print(clone_str)
  
  # chrs <- c(as.character(1:22), "X", "Y")
  chrs <- c(as.character(1:22), "X")
  df_track <- annotables::grch37 %>% 
    dplyr::select(gene_symbol = symbol, chr, start, end) %>% 
    inner_join(deg_df) %>% 
    dplyr::filter(chr %in% chrs) %>% 
    dplyr::mutate(chr = factor(chr, levels = chrs),
                  position = (start + end) / 2)
  
  # MA: sometimes a genes appears more than once, keep only the unique ones
  df_track <- distinct(df_track)
  # cols <- rev(brewer.pal(n = 3, name = "RdBu") )
  # cols <- cols[1:length(unique(deg_df$is_genome))]
  cols <- rev(brewer.pal(n = 6, name = "RdBu"))
  # library(scales)

  View(head(df_track))
  colnames(df_track)
  dim(df_track)
  # df_track[nrow(df_track)+1,] <- c('NA','14','0','0',rep('NA',17))
  track_plot <- ggplot(df_track, aes(x = position, y = mlog10FDR)) +  #-log10(p_val_adj)
    geom_point(aes(colour = logFC), size = 1.3) +   # MA: size was 1
    facet_grid(~ chr, scales = "free_x", space='free',
               # strip.position = "bottom",
               switch = "x",
               # nrow = 1
    ) +
    scale_colour_gradientn(name = "logFC",
                           colours = cols,
                           values = rescale(c(-3, -1.5, -0.25, 0, 0.25, 1.5, 3)),
                           limits = c(-4.25, 4.25)) +
    theme(#strip.background = element_blank(),
      strip.background = element_rect(fill = 'white', colour = 'white'),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=10),   
      legend.text = element_text(size=10),
      legend.title = element_text(size=10),
      strip.placement = "outside",
      legend.position = "top",
      panel.spacing = unit(0.1, 'cm'),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # panel.background = element_blank(), 
      axis.line = element_line(colour = "black")) +   ## MA: was 0.2
    theme(text = element_text(size = 18)) +
    # theme_classic() +
    labs(x = "Chromosome",y = "-Log10(FDR)",
         subtitle = clone_str) +
    scale_x_continuous(expand = c(0,0))
  
  
  
  df_annot <- top_n(df_track, n_genes_to_annotate, abs(logFC))
  
  
  print("The labeled genes")
  print(dim(df_annot))
  
  if(!is.null(additional_genes)) {
    df_annot <- df_annot  %>% 
      bind_rows(dplyr::filter(df_track, gene_symbol %in% additional_genes))
  }
  
  ## this is where the genes are added
  track_plot2 <- track_plot +
    geom_text_repel(data = df_annot, aes(label = gene_symbol), size = 2.4)
  
  # base_name <- paste(obs_clones, collapse = "_vs_")
  # png(paste0(save_dir,"DE_chr_",base_name,".png"), height = 2*450, width=2*1200,res = 2*72)
  # print(track_plot2)
  # dev.off()
  
  return(track_plot2)
}

# df_cnv has column names: ensembl_gene_id
plot_CNV <- function(df_cnv_fn, clones, meta_genes=NULL){
  df_cnv <- data.table::fread(df_cnv_fn, check.names = F, stringsAsFactors = F) %>% as.data.frame()
  df_cnv <- df_cnv %>%
    dplyr::rename(ensembl_gene_id=V1)
  
  
  print(dim(df_cnv))
  # Select only genes that contain variance across all clones --> make plot clearer
  if(is.null(meta_genes)){
    var_genes <- df_cnv %>%
      tibble::column_to_rownames('ensembl_gene_id') %>%
      dplyr::mutate(var_gene = rowVars(as.matrix(.), na.rm=T)) %>%
      dplyr::pull(var_gene)
    
    df_cnv <- df_cnv[var_genes>0,]
  }else{
    df_cnv <- df_cnv %>%
      dplyr::filter(ensembl_gene_id %in% meta_genes$ensembl_gene_id)
  }
  
        
  df_cnv <- df_cnv %>%
    # dplyr::rename(ensembl_gene_id=V1)%>%
    dplyr::select(all_of(c('ensembl_gene_id',clones)))%>%
    tidyr::pivot_longer(!ensembl_gene_id, values_to = 'cnv', names_to = 'clone')%>%
    as.data.frame()
  # df_cnv <- df_cnv %>% 
  #   dplyr::filter(clone %in% clones)  #, !chr %in% c("X","Y")
  print(dim(df_cnv))
  # sum(is.na(df_cnv$cnv))
  df_cnv <- df_cnv[!is.na(df_cnv$cnv),]
  # df_cnv1 <- df_cnv
  annots <- annotables::grch38
  annots <- dplyr::select(annots, ensembl_gene_id = ensgene, chr, start, end) 
  ## MA: 24 Jan 2021: remove duplicated rows, these would result in spaces in the CN track
  annots <- annots[!duplicated(annots),]
  # backup_cnv <- df_cnv
  # df_cnv <- backup_cnv
  genes1 <- df_cnv[df_cnv$clone==clones[1],'ensembl_gene_id']
  genes2 <- df_cnv[df_cnv$clone==clones[2],'ensembl_gene_id']
  df_cnv <- df_cnv[df_cnv$ensembl_gene_id %in% intersect(genes1,genes2),]
  print(dim(df_cnv))
  
  df_cnv <- inner_join(df_cnv, annots)
  print(dim(df_cnv))
  chr_levels <- c(as.character(1:22), "X")
  df_cnv <- df_cnv %>% dplyr::filter(chr %in% chr_levels)
  # t <- df_cnv %>% dplyr::filter(chr ==17)
  # dim(t)
  # dim(df_track)
  # t$ensembl_gene_id
  # t <- t %>% dplyr::filter(ensembl_gene_id %in% df_track$ensembl_gene_id)
  # summary(as.factor(df_track$chr))
  
  df_cnv <- dplyr::select(df_cnv, cnv, clone, ensembl_gene_id, chr, start, end) %>%
    group_by(chr) %>%
    # dplyr::mutate(start_order = rank((start+end)/2)) %>%
    dplyr::mutate(start_order = rank(start)) %>%
    ungroup() 
    # %>%
    # inner_join(df_cnv)
  # df_cnv$position <- (df_cnv$start+df_cnv$end)/2
  # genes_use <- union(var_genes, intersect(rownames(cnv), df_cnv$ensembl_gene_id))
  # length(genes_use)
  # df_cnv <- df_cnv %>%
  #      dplyr::filter(ensembl_gene_id %in% genes_use)
  
  
  # length(unique(df_cnv$ensembl_gene_id))
  # df_cnv1 <- df_cnv %>%
  #   dplyr::group_by(ensembl_gene_id) %>%
  #   dplyr::summarise(is_variant=var(cnv)>0) %>%
  #   dplyr::filter(is_variant==T)
  # 
  # df_cnv <- df_cnv %>%
  #   dplyr::filter(ensembl_gene_id %in% unique(df_cnv1$ensembl_gene_id))
  # length(unique(df_cnv$ensembl_gene_id))
  
  # print(summary(as.factor(df_cnv$clone)))
  # df_cnv$start[1:3]
  
  # df_cnv <- dplyr::select(df_cnv, ensembl_gene_id, chr, start, end) %>% 
  #   group_by(chr) %>% 
  #   dplyr::mutate(start_order = rank((start + end) / 2)) %>% # #
  #   ungroup() %>% 
  #   inner_join(df_cnv)
  
  
  
  # chr_levels <- c(as.character(1:23), "X", "Y")
  df_cnv$chr <- factor(df_cnv$chr, levels = chr_levels)
  # summary(as.factor(df_cnv$chr))
  # dim(df_cnv)
  df_cnv <- tidyr::drop_na(df_cnv)
  # df_cnv1 <- df_cnv
  # MA: 11 Apr 2020, setting the copy number colors that were used in the heatmap
  # TO COME BACK
  #print("data frame")
  #cnv_colors <- c('#4880B8', '#A7C9DF','#CCCCCC','#F5CE93','#ED9364','#D2553E','#A42116','#8B1A43','#CB3576','#D06CAD','#C196C4','#D0BAD8')
  #cnv_cols <- data.frame(cn=0:11,color<-cnv_cols)
  #colnames(cnv_cols) <- c("cn","color")
  #levels(cnv_cols$color) <- cnv_colors
  
  cnv_colors <- c('#4880B8', '#A7C9DF','#CCCCCC','#F5CE93','#ED9364','#D2553E','#A42116','#8B1A43','#CB3576','#D06CAD','#C196C4','#D0BAD8')
  
  cnv_cols <- c('0'='#4880B8', '1'='#A7C9DF','2'='#CCCCCC','3'='#F5CE93','4'='#ED9364',
                '5'='#D2553E','6'='#A42116','7'='#8B1A43','8'='#CB3576','9'='#D06CAD',
                '10'='#C196C4','11'='#D0BAD8')
  #levels(cnv_cols) <- 0:11
  #  cnv_cols <- c("0" = "#2166ac",
  #                "1" = "#92c5de", 
  #                "2" = "grey80", 
  #                "3" = "#f4a582", 
  #                "4" = "#d6604d",
  #                "5" = "#b2182b",
  #                "6+" = "#67001f")
  
  df_cnv$cnv <- as.character(round(df_cnv$cnv))
  levels(df_cnv$cnv) <- 0:(length(cnv_cols)-1)
  
  #print(levels(df_cnv$cnv))
  
  # MA: removing the 6+ restriction (clonealign still has that restriction though)
  # df_cnv$cnv[round(df_cnv$median_cnmode) >= 6] <- "6+"
  #!chr %in% c("X","Y")
  
  # Just in case of SA609
  df_cnv$clone <- ifelse(df_cnv$clone=='R','A',df_cnv$clone)
  if(clones[1]=='R'){
    clones[1] <- 'A'
  }
  
  df_cnv$clone <- factor(df_cnv$clone, levels = c(clones[2],clones[1]))
  print(summary(df_cnv$clone))
  # saveRDS(df_cnv, paste0(save_fig_dir,'df_cnv_1035_test.rds'))
  cnv_plot <- df_cnv %>% 
    ggplot(aes(x = start_order, y = clone, fill = cnv)) + #, fill = cnv
    # geom_tile(aes(fill = cnv)) + #
    geom_raster() +
    # geom_tile() + 
    facet_wrap(~ chr, scales = "free_x", nrow = 1, switch = "x") +
    # facet_grid(~ chr, scales = "free_x", switch = "x", space='free') +
    # theme(legend.position = "bottom", axis.text.x = element_blank()) +
    scale_y_discrete(expand = c(0, 0)) +
    # scale_fill_manual(values=cnv_colors, name = "Copy number", guide = 'legend',labels = 0:(length(cnv_colors)-1),drop=FALSE)  +
    #          theme(legend.position = "bottom") +
    scale_fill_manual(values = cnv_cols, name = "Copy number ", breaks = names(cnv_cols)) +  #
    labs(x = "Chromosome", y = "Clone")
  
    my_font <- "Helvetica"
  
    lg_pos <- "bottom"
    # thesis_theme <- ggplot2::theme(
    #   text = element_text(color="black",size = 8, hjust = 0.5, family=my_font),
    #   axis.title.x = element_text(color="black",size=8, hjust = 0.5, family=my_font),
    #   axis.title.y = element_text(color="black",size=8, hjust = 0.5, family=my_font),
    #   # axis.text.x = element_text(color="black",size=7, hjust = 0.5, family=my_font, angle = 90),
    #   axis.text.x = element_blank(),
    #   axis.ticks.x = element_blank(),
    #   strip.placement = "outside",
    #   axis.line = element_line(colour = "black"),
    #   axis.text.y = element_text(color="black",size=7, hjust = 0.5, family=my_font),
    #   plot.title = element_text(color="black",size=10, face="bold", hjust=0, family=my_font),
    #   legend.title=element_text(color="black",size=7, hjust = 0.5, family=my_font),
    #   legend.text=element_text(color="black",size=7, hjust = 0.5, family=my_font),
    #   strip.text.x = element_text(color="black",size=9, family=my_font),
    #   strip.text.y = element_text(color="black",size=9, family=my_font),
    #   # legend.spacing.x = unit(0.1, 'mm'),
    #   # legend.spacing.y = unit(0.1, 'mm'),
    #   # legend.key.height=unit(1,"line"),
    #   legend.position = lg_pos,
    #   panel.background = element_blank(), 
    #   panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    # )
    # cnv_theme <- theme(strip.background = element_rect(fill = 'white', colour = 'white'),
    #                    axis.text.x = element_blank(),
    #                    axis.ticks.x = element_blank(),
    #                    axis.text.y = element_text(size=12, colour = "black", hjust = 0.5),
    #                    axis.title.y = element_text(size=12, colour = "black", hjust = 0.5),
    #                    axis.line = element_line(colour = "black"),
    #                    strip.placement = "outside",
    #                    legend.position = "bottom",
    #                    legend.text = element_text(size=9),
    #                    legend.title = element_text(size=9),
    #                    legend.key.size=unit(0.3,"cm"),
    #                    panel.grid.major = element_blank(),
    #                    panel.grid.minor = element_blank(),
    #                    # panel.background = element_rect(fill = "#F8F8F8", colour = NA),
    #                    panel.spacing = unit(c(0.1), 'cm'))
    cnv_plot <- cnv_plot + theme(strip.background = element_rect(fill = 'white', colour = 'white'),
                       text = element_text(color="black",size = 8, hjust = 0.5, family=my_font),
                      axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.text.y = element_text(color="black",size=7, hjust = 0.5, family=my_font),
                      axis.title.y = element_text(color="black",size=8, hjust = 0.5, family=my_font),
                      axis.line = element_line(colour = "black"),
                      strip.placement = "outside",
                      legend.position = "bottom",
                      legend.text=element_text(color="black",size=7, hjust = 0.5, family=my_font),
                      legend.title=element_text(color="black",size=7, hjust = 0.5, family=my_font),
                      legend.key.size=unit(0.3,"cm"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      # panel.background = element_rect(fill = "#F8F8F8", colour = NA),
                      panel.spacing = unit(c(0.1), 'cm'),
                      legend.margin=margin(0,0,0,0),
                      legend.box.margin=margin(-2,-2,-2,-2))    # MA: was 0.2

    cnv_plot <- cnv_plot + guides(fill = guide_legend(nrow = 1, override.aes = list(size=0.1))) +  #, override.aes = list(size=1.1)
                           scale_x_continuous(expand = c(0,0))
  
  # cnv_plot
  lg <- cowplot::get_legend(cnv_plot)
  plg <- cowplot::ggdraw() + cowplot::draw_plot(lg)
    
  results <- list(df_cnv=df_cnv, plg=plg, cnv_plot=cnv_plot)
  return(results)
}



viz_trackplot <- function(save_dir){
  main_plot <- cowplot::plot_grid(
    track_plot2 + theme(axis.title.x = element_blank(),
                        strip.text.x = element_blank()),
    cnv_plot,
    ncol = 1,
    rel_heights = c(3,1),
    align = 'v'
  )
  save_dir <- input_dir
  png(paste(save_dir,"track_plots.png",sep="/"), height = 2*500, width=2*1200, res = 2*72)
  print(main_plot)
  dev.off()
  
  
  return(main_plot)
}




