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

plot_gex_cnv_v2 <- function(deg_df,
                         df_cnv,
                         save_dir,
                         plt_subtitle='',
                         clone_str=NULL,
                         additional_genes_fn=NULL,
                         n_genes_to_annotate=40,
                         plttitle=NULL,
                         gene_type='in-cis',
                         deg_fn=NULL) {
  # print(deg_fn)
  # plttitle <- clone_str
  if(is.null(plttitle)){
    tl <- element_blank()
  }else{
    tl <- element_text(size=15, colour = "black", face='bold')  # 
  }
  if(is.null(clone_str)){
    # clone_str <- paste(clones, collapse = " vs ")
    clone_str <- ''
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
  if(is.null(deg_df)){
    deg_df <- read.csv(deg_fn, check.names=F, stringsAsFactors=F) #, row.names = 1  
  }
  print(dim(deg_df))
  deg_df$chr <- NULL # bulk data, duplicated chr column
  # deg_df <- deg_df %>%
  #   dplyr::filter(abs(logFC)>0.25 & FDR<0.01 & PValue<0.05)
  # print(dim(deg_df))
  
  print(dim(deg_df))
  if('ens_gene_id' %in% colnames(deg_df)){
    deg_df <- deg_df %>%
      dplyr::rename(ensembl_gene_id=ens_gene_id)
  }
  chrs <- c(as.character(1:22), "X")
  # annots <- dplyr::select(annots, ensembl_gene_id = ensgene, chr, start, end) 
  df_track <- annotables::grch37 %>%
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
  
  
  cols <- rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu"))
  # cols <- c("#021621", "#92C5DE", "#F7F7F7", "#F4A582", "#610614")
  #blue: #053885, red: #8f0319
  # df_track$mlog10FDR <- -log10(df_track$FDR)
  # df_track$mlog10FDR <- sapply(df_track$mlog10FDR, function(x) replace(x, is.infinite(x), 270))
  # df_track$mlog10FDR <- sapply(df_track$mlog10FDR, function(x) replace(x, x>270, 270))
  max_logFC <- 3
  df_track$logFC <- sapply(df_track$logFC, function(x) replace(x, x > max_logFC, max_logFC))
  df_track$logFC <- sapply(df_track$logFC, function(x) replace(x, x < (-max_logFC), -max_logFC))
  if(min(df_track$logFC)<(-3)){
    ymin <- -3
  }else{
    ymin <- min(df_track$logFC)
  }
  if(max(df_track$logFC)<3){
    ymax <- 3
  }else{
    ymax <- max(df_track$logFC)
  }
  df_cnv$chr <- as.character(df_cnv$chr)
  blank_data <- NULL
  for (chr in unique(df_cnv$chr)) {
    # print(min(df_cnv[df_cnv$chr==chr,'start_order']))
    # print(max(df_cnv[df_cnv$chr==chr,'start_order']))
    blank_data <- rbind(blank_data, data.frame(chr=chr, start_order=min(df_cnv[df_cnv$chr==chr,'start_order']), logFC=0))
    blank_data <- rbind(blank_data, data.frame(chr=chr, start_order=max(df_cnv[df_cnv$chr==chr,'start_order']), logFC=0))
  }
  print("Blank data")
  # print(blank_data)
  
  # df_track1 <- df_track
  # df_track <- df_track1
  df_track$start_order <- 0.0
  genes_use <- intersect(df_track$ensembl_gene_id, df_cnv$ensembl_gene_id)
  length(genes_use)
  for (geneid in genes_use) {
    t <- df_cnv[df_cnv$ensembl_gene_id==geneid,]$start_order
    df_track[df_track$ensembl_gene_id==geneid,'start_order'] <- as.numeric(t[1])
  }
  # df_track <- drop_na(df_track)
  dim(df_track)
  df_track$chr <- factor(df_track$chr, levels = chrs)
  # track_plot <- ggplot(df_track, aes(x = position, y = logFC)) +  #-log10(FDR)
  #   geom_point(aes(colour = logFC), size = 1.5) +   # MA: size was 1 , shape=gene_type
  track_plot <- ggplot(df_track, aes(x = start_order, y = logFC)) +
    geom_point(aes(colour = logFC), size = 1.3) +   # MA: size was 1
    geom_blank(data = blank_data, aes(x = start_order, y = logFC)) +
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
    scale_colour_gradientn(name = '  log2FC ',  #paste0(clone_str,'  log2FC ')
                           colours = cols, 
                           values = rescale(c(-3, -1, 0, 1, 3)), #rescale(c(-2, -0.3, 0, 0.3, 2)
                           limits = c(-3, 3)) +
    ylim(ymin, ymax)
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
      legend.key.height=unit(1,"line"),
      legend.position = lg_pos,
      panel.background = element_blank(), 
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    )
   track_plot <- track_plot + thesis_theme
   track_plot <- track_plot + labs(x = "Chromosome", y="log2 FC", title = plttitle, subtitle = plt_subtitle) +  #subtitle = sample_name
                              scale_x_continuous(expand = c(0,0))
  
  # track_plot
  # df_annot <- top_n(df_track, n_genes_to_annotate, -log10(FDR))
  # df_annot <- top_n(df_track, n_genes_to_annotate, logFC) 
  # df_track <- df_track[order(df_track$logFC, df_track$mlog10FDR, decreasing = T),] 
  df_track <- df_track[order(abs(df_track$logFC), decreasing = T),] 
  # 
  df_annot <- df_track[1:n_genes_to_annotate,]
  print(df_annot$gene_symbol)
  # dim(df_annot)
  # unique(df_track$chr)
  # dim(df_track)
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
    geom_text_repel(data = df_annot, aes(label = gene_symbol), size = 3.1,# was 3.3
                    max.overlaps = Inf,
                    min.segment.length = 0,
                    ylim = c(-Inf, Inf), xlim= c(-Inf, Inf),
                    color='black', segment.alpha = 0.2)
                    #nudge_x = 0.1, direction = "x", hjust = "right",
                    #segment.ncp = 3)
                    # segment.curvature = -0.1,
                    # )
                    #family = my_font)
  
  
  return(track_plot2)
}




# df_cnv has column names: ensembl_gene_id
plot_CNV <- function(df_cnv, sample_ids, mainsites, sample_tag='SA919X7XB0'){ #, df_cnv_fn=NULL
  # df_cnv <- read.csv(df_cnv_fn, check.names = F, stringsAsFactors = F)
  chrs <- c(as.character(1:22), "X")
  df_cnv <- df_cnv %>%
    tidyr::pivot_longer(!ensembl_gene_id, names_to = "clone", values_to = "cnv")
  print(unique(df_cnv$clone))
  # sample_ids <- c(as.character(s1),as.character(s2))
  sample_ids <- gsub(sample_tag,'',sample_ids)
  sample_ids <- paste0(mainsites,'_',sample_ids)
  df_cnv$clone <- ifelse(df_cnv$clone=='s1',sample_ids[1],sample_ids[2])
  # df_cnv$clone <- gsub('SA919X7XB0','',df_cnv$clone)
  # df_cnv$clone <- ifelse(df_cnv$clone=='5604',paste0('Met_',df_cnv$clone),paste0('Pri_',df_cnv$clone))
  df_cnv$clone <- factor(df_cnv$clone, levels = c(sample_ids[2],sample_ids[1]))
  
  # df_cnv <- df_cnv %>% 
  #   dplyr::filter(clone %in% clones)  #, !chr %in% c("X","Y")
  dim(df_cnv)
  annots <- annotables::grch37 %>%
     dplyr::select(ensembl_gene_id = ensgene, chr, start, end) %>%
     dplyr::filter(chr %in% chrs)
  ## MA: 24 Jan 2021: remove duplicated rows, these would result in spaces in the CN track
  annots <- annots[!duplicated(annots$ensembl_gene_id),]
  # backup_cnv <- df_cnv
  # df_cnv <- backup_cnv
  # genes1 <- df_cnv[df_cnv$clone==clones[1],'ensembl_gene_id']
  # genes2 <- df_cnv[df_cnv$clone==clones[2],'ensembl_gene_id']
  # common <- intersect(genes1,genes2)
  # length(common)
  # df_cnv <- df_cnv[df_cnv$ensembl_gene_id %in% common,]
  # print(dim(df_cnv))
  
  df_cnv <- inner_join(df_cnv, annots, by=c('ensembl_gene_id'))
  print(dim(df_cnv))
  
  
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
  
  chr_levels <- c(as.character(1:22), "X")
  df_cnv <- df_cnv %>% dplyr::filter(chr %in% chr_levels)
  
  
  
  
  # chr_levels <- c(as.character(1:23), "X", "Y")
  df_cnv$chr <- factor(df_cnv$chr, levels = chr_levels)
  # summary(as.factor(df_cnv$chr))
  df_cnv <- tidyr::drop_na(df_cnv)
  dim(df_cnv)
  
  # cnv_colors <- c('#4880B8', '#A7C9DF','#CCCCCC','#F5CE93','#ED9364','#D2553E','#A42116','#8B1A43','#CB3576','#D06CAD','#C196C4','#D0BAD8')
  
  cnv_cols <- c('0'='#4880B8', '1'='#A7C9DF','2'='#CCCCCC','3'='#F5CE93','4'='#ED9364',
                '5'='#D2553E','6'='#A42116','7'='#8B1A43','8'='#CB3576','9'='#D06CAD',
                '10'='#C196C4','11'='#D0BAD8')

  df_cnv$cnv <- as.character(round(df_cnv$cnv))
  levels(df_cnv$cnv) <- 0:(length(cnv_cols)-1)
  
  #print(levels(df_cnv$cnv))
  
  # MA: removing the 6+ restriction (clonealign still has that restriction though)
  # df_cnv$cnv[round(df_cnv$median_cnmode) >= 6] <- "6+"
  #!chr %in% c("X","Y")
  # df_cnv$clone <- ifelse(df_cnv$clone=='R','A',df_cnv$clone)
  # if(clones[1]=='R'){
  #   clones[1] <- 'A'
  # }
    # df_cnv$clone <- factor(df_cnv$clone, levels = c(clones[2],clones[1]))
  # 
  # df_cnv$clone <- factor(df_cnv$clone, levels = c(clones[2],clones[1]))
  # print(summary(as.factor(df_cnv$clone)))
  
  # saveRDS(df_cnv, paste0(save_fig_dir,'df_cnv_1035_test.rds'))
  cnv_plot <- df_cnv %>% 
    ggplot(aes(x = start_order, y = clone, fill = cnv)) + #, fill = cnv
    # geom_tile(aes(fill = cnv)) + #
    geom_raster() +
    # geom_tile() + 
    facet_wrap(~ factor(chr, levels = chrs), scales = "free_x", nrow = 1, switch = "x", drop=F) +
    # facet_grid(~ chr, scales = "free_x", switch = "x", space='free') +
    # theme(legend.position = "bottom", axis.text.x = element_blank()) +
    scale_y_discrete(expand = c(0, 0)) +
    # scale_fill_manual(values=cnv_colors, name = "Copy number", guide = 'legend',labels = 0:(length(cnv_colors)-1),drop=FALSE)  +
    #          theme(legend.position = "bottom") +
    scale_fill_manual(values = cnv_cols, name = "Copy number ", breaks = names(cnv_cols)) +  #
    labs(x = "Chromosome", y = "Sample") #+ 
    # ylim(ymin, ymax)
  
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
                      panel.spacing = unit(c(0.1), 'cm'))    # MA: was 0.2

  cnv_plot <- cnv_plot + guides(fill = guide_legend(nrow = 1, override.aes = list(size=0.1))) +  #, override.aes = list(size=1.1)
                           scale_x_continuous(expand = c(0,0))
  
  # cnv_plot
  results <- list(df_cnv=df_cnv, cnv_plot=cnv_plot)
  return(results)
}



viz_trackplot <- function(res, save_dir, comps, sample_tag='SA919X7XB0'){
  obs_cp <- res$comparison
  rownames(comps) <- comps$fn
  sample_ids <- c(res$s2, res$s1)
  mainsites <- c(comps[obs_cp,'mainsite_s1'],comps[obs_cp,'mainsite_s2'])
  mainsites <- gsub('Metastasis','Met',mainsites)
  mainsites <- gsub('Primary','Pri',mainsites)
  # comps <- data.table::fread(paste0(bulk_dir,'comparison_bulk_SA919.csv'))%>% as.data.frame()
  
  res_cnv <- plot_CNV(res$cnv_total, sample_ids, mainsites, sample_tag)
  desc <- summary(as.factor(res$cnv$positive_change))
  linear_trend <- ''
  for(x in names(desc)){
    linear_trend <- paste0(linear_trend, x,': ',round(desc[x]/dim(res$cnv)[1]*100,1),'%; ')  
  }
  
  # save_dir <- paste0(bulk_dir,obs_cp,'/')
  # plttitle <- paste0(mainsites[1],':',gsub(sample_tag,'',res$s1),' vs. ',mainsites[2],':',gsub(sample_tag,'',res$s2))
  plttitle <- 'Metastasis vs Primary'
  pct_change <- round(comps[obs_cp,'nb_signf']/comps[obs_cp,'nb_cnchange']*100,1)
  plt_subtitle <- paste0('Proportion of CN change coupled with gene exp: ',pct_change,'%\n',linear_trend)
  cis_de <- res$cnv
  class(cis_de)
  sum(is.na(res_cnv$df_cnv))
  print(summary(as.factor(cis_de$chr)))
  # cis_de <- cis_de
  #   dplyr::filter()
  
  track_plot2 <- plot_gex_cnv_v2(res$cnv,
                 res_cnv$df_cnv,
                 save_dir,
                 plt_subtitle,
                 clone_str=NULL,
                 additional_genes_fn = NULL,
                 n_genes_to_annotate = 20,
                 plttitle,
                 gene_type='in-cis',
                 deg_fn=NULL)
  # track_plot2  
    
  main_plot <- cowplot::plot_grid(
    track_plot2 + theme(axis.title.x = element_blank(),
                        strip.text.x = element_blank()),
    res_cnv$cnv_plot,
    ncol = 1,
    rel_heights = c(4,1),
    align = 'v'
  )

  # png(paste0(save_dir,obs_cp,"_track_plots.png"), height = 2*700, width=2*1300, res = 2*72)
  # print(main_plot)
  # dev.off()
  
  ggsave(paste0(save_dir,"track_plots.png"),
         plot = main_plot,
         height = 7,
         width = 13,
         # useDingbats=F,
         type = "cairo-png",
         dpi=150
  )
  # return(main_plot)
}




